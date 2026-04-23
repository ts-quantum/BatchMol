import numpy as np
import os, sys
base_dir = os.path.dirname(sys.executable)
if base_dir not in sys.path:
    sys.path.insert(0, base_dir)

# set path for Matplotlib config (e.g. ~/.matplotlib_cache)
mpl_cache_dir = os.path.expanduser('~/.matplotlib_cache')
if not os.path.exists(mpl_cache_dir):
    os.makedirs(mpl_cache_dir)

# set env var Mpl_CONFIGDIR 
os.environ['MPLCONFIGDIR'] = mpl_cache_dir

import matplotlib.pyplot as plt
from matplotlib import colormaps
from contextlib import redirect_stderr
from tqdm import tqdm # progress bar!!
import pyscf.tools.molden as molden_tools
from cclib.io import ccread
import pyvista as pv
import psutil 
import platform,subprocess
from pyscf import scf, lib, data, dft, gto 
from collections import defaultdict
import click # input during execution

"""
BUILD - macOS
python -m nuitka --standalone --macos-create-app-bundle --macos-app-name="BatchMol" --enable-plugin=numpy --enable-plugin=matplotlib --enable-plugin=anti-bloat --nofollow-import-to=pyscf --nofollow-import-to=vtkmodules --no-deployment-flag=excluded-module-usage --no-deployment-flag=self-execution --jobs=8 --output-dir=dist --remove-output batch.py

"""

class MoleculeData: #always call arguments by name! 
                    #no special order "*"
                    #e.g. obj=MoleculeData(atoms=atoms, type="molden")
    def __init__(self, *, name=None, atom_points=None, atom_types=None,
                 mol=None, mo_energy=None, mo_coeff=None, mo_occ=None,
                 orb_labels=None, spin=None, esp_mesh=None
                 ):
        self.name=name
        self.atom_points=atom_points
        self.atom_types=atom_types
        self.mol=mol
        self.mo_energy=mo_energy
        self.mo_coeff=mo_coeff
        self.mo_occ=mo_occ
        self.orb_labels=orb_labels
        self.spin=spin
        self.esp_mesh=esp_mesh

    @classmethod
    def from_molden(cls, filepath):
        short_name = os.path.basename(filepath)
        with open(os.devnull, 'w') as devnull:
                with redirect_stderr(devnull):
                    mol, mo_energy, mo_coeff, mo_occ, orb_labels, spin = molden_tools.load(filepath)
        mol.cart = False #####
        mol.build()
        atom_points = mol.atom_coords() * 0.529177 
        atom_types = [data.elements.charge(mol.atom_symbol(i))
                        for i in range(mol.natm)]
        return cls(
            name=short_name, 
            atom_points=atom_points,
            atom_types=atom_types,
            mol=mol,
            mo_energy=mo_energy,
            mo_coeff=mo_coeff,
            mo_occ=mo_occ,
            orb_labels=orb_labels,
            spin=spin,
            esp_mesh=None
        )
    @classmethod
    def fix_fchk_format(cls,filepath):
        import re
        import io    
        with open(filepath, 'r') as f:
            lines = f.readlines()
        fixed_lines = []
        for line in lines:
            # check and repair all lines containing data (if necessary)
            if line.startswith(" "):
                # look for number followed by sign (no e/E)
                # e.g.: 0.123-4.567 -> 0.123 -4.567
                line = re.sub(r'(\d)([+-]\d\.)', r'\1 \2', line)
                # Special case for 'lost' spaces with extreme exponents.
                # e.g.: 5.78e-02-4.48e-149 -> 5.78e-02 -4.48e-149
                line = re.sub(r'(e[+-]\d{2,3})([+-]\d\.)', r'\1 \2', line)
                
            fixed_lines.append(line)
            
        return io.StringIO("".join(fixed_lines))

    @classmethod
    def from_fchk(cls, filepath):
        
        # Apply fix
        fixed_file_stream = MoleculeData.fix_fchk_format(filepath)
        # 1. load data with cclib 
        data = ccread(fixed_file_stream)
        
        # 2. manually build PySCF molecule
        mol = gto.Mole()
        # Atoms from cclib (atomnos = atomic number, atomcoords in Angström)
        mol.atom = [[data.atomnos[i], data.atomcoords[-1][i]] for i in range(data.natom)]
        
        mol.basis = data.metadata.get('basis_set')

        mol.charge = data.charge
        mol.spin = data.mult - 1
        
        mol.cart = True 
        mol.build()

        n_ao_fchk = data.mocoeffs[0].shape[1]
        
        if mol.nao != n_ao_fchk:
            mol.cart = False
            mol.build()
            if mol.nao != n_ao_fchk:
                print(f"KRITISCH: Basis-Mismatch! PySCF AOs: {mol.nao}, FCHK AOs: {n_ao_fchk}")

        # 4. Data Extraction (Alpha/Beta Handling)
        is_uhf = len(data.homos) == 2
        if is_uhf:
            mo_coeff = [m.T for m in data.mocoeffs] # Transpose for (AO, MO)
            mo_energy = data.moenergies
            # Generate Occupation
            mo_occ = [np.zeros(len(e)) for e in mo_energy]
            mo_occ[0][:data.homos[0]+1] = 1.0
            mo_occ[1][:data.homos[1]+1] = 1.0
            orb_labels = [[f"Alpha MO {i}" for i in range(len(mo_energy[0]))],
                        [f"Beta MO {i}" for i in range(len(mo_energy[1]))]]
        else:
            mo_coeff = data.mocoeffs[0].T
            mo_energy = data.moenergies[0]
            mo_occ = np.zeros(len(mo_energy))
            mo_occ[:data.homos[0]+1] = 2.0
            orb_labels = [f"MO {i}" for i in range(len(mo_energy))]

        return cls(
            name=os.path.basename(filepath),
            atom_points=data.atomcoords[-1],
            atom_types=data.atomnos.tolist(),
            mol=mol,
            mo_energy=mo_energy,
            mo_coeff=mo_coeff,
            mo_occ=mo_occ,
            orb_labels=orb_labels,
            spin=mol.spin
        )

################### Plotting #####################################
def get_radius_by_group(atomic_number):
    # Definition of periods (Start, End): Radius
    groups = {
        (1, 1): 0.2,    # Hydrogen
        (2, 2): 0.2,    # Helium
        (3, 10): 0.35,   # 2. Period (Li to Ne)
        (11, 18): 0.45,  # 3. Period (Na to Ar)
        (19, 36): 0.55,  # 4. Period (K to Kr)
        (37, 54): 0.65,  # 5. Period (Rb to Xe)
        (55, 86): 0.75,  # 6. Period (Cs-Rn, incl. Pt, Au)
    }
    
    for (start, end), radius in groups.items():
        if start <= atomic_number <= end:
            return radius
    return 0.3  # Default

def draw_mol(atom_points, atom_types, cpk_colors):
    all_parts = []
    # create PolyData-Object from all points
    atoms_poly = pv.PolyData(atom_points)
    # add atom-tpye as scalar for coloring
    atoms_poly.point_data["colors"] = atom_types
    #sphere as template
    #sphere_source = pv.Sphere(radius=0.3, theta_resolution=20, phi_resolution=20)
    # color mapping:
    # glyph object contains original atom IDs, Lookup Table (LUT) can be used
    # loop over type  
    u_types = np.unique(atom_types)
    for atom_type in u_types:
        color = cpk_colors[atom_type]
        mask = atom_types == atom_type
        if np.any(mask):
            sub_atoms = atoms_poly.extract_points(mask)
            r=get_radius_by_group(atom_type)
            sphere_source = pv.Sphere(radius=r, theta_resolution=20, phi_resolution=20)   
            glyphs = sub_atoms.glyph(geom=sphere_source, scale=False, orient=False)
            
            rgb = (np.array(pv.Color(color).float_rgb) * 255).astype(np.uint8)
            colors_array = np.tile(rgb, (glyphs.n_points, 1))
            glyphs.point_data["RGB"] = colors_array
            all_parts.append(glyphs)
    
    # --- Bonds as single net
    lines = []
    for i in range(len(atom_points)):
        type_i = int(atom_types[i])
        rad_i = cov_radii.get(type_i, default_radius)
        for j in range(i + 1, len(atom_points)):
            type_j = int(atom_types[j])
            rad_j = cov_radii.get(type_j, default_radius) 
            dist = np.linalg.norm(atom_points[i] - atom_points[j])
            bd_threshold = rad_i + rad_j + 0.6
            if 0.6 < dist < bd_threshold:
                # only indices of linked points are saved
                lines.append([2, i, j]) # 2: line conects two points      
    tubes = None
    if lines:
        # Create PolyData-Object for lines
        bonds_poly = pv.PolyData(atom_points)
        bonds_poly.lines = np.hstack(lines)
        # convert lines into tubes
        tubes = bonds_poly.tube(radius=0.06)
        # bond color
        bond_rgb = (np.array(pv.Color("lightgray").float_rgb) * 255).astype(np.uint8)
        tubes.point_data["RGB"] = np.tile(bond_rgb, (tubes.n_points, 1))
        all_parts.append(tubes)
    # merge all meshes
    if not all_parts:
        return None 
    combined = all_parts[0].merge(all_parts[1:])
    return combined

def draw_orb_molden(data_obj, orbital_index=0, spin_idx=0, iso_level=0.02,
                    nx = 50, ny = 50, nz = 50, padding = 3.0):
    visual_objects = [] # list object for mesh and color
    # extract coordinates
    # PySCF uses Bohr, Angstrom is needed for plotting
    atom_coords_angstrom = data_obj.mol.atom_coords() * 0.529177 # conversion Bohr -> Angström

    # set Grid boundaries (in Angström)
    padding = padding
    xmin, ymin, zmin = np.min(atom_coords_angstrom, axis=0) - padding
    xmax, ymax, zmax = np.max(atom_coords_angstrom, axis=0) + padding

    # calculate Grid points (resolution nxxnyxnz)
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    z = np.linspace(zmin, zmax, nz)

    # PySCF eval_gto uses coordinates in Bohr
    coords_angstrom = lib.cartesian_prod([x, y, z])
    coords_bohr = coords_angstrom / 0.529177

    # GTO-Werte at points (AO = Atomic Orbitals)
    # "GTOval" selects automatically "spherical" or "cartesian" depending on mol object
    ao_values = data_obj.mol.eval_gto("GTOval", coords_bohr)

     # Select the correct MO coefficients (handle RKS vs UKS/UHF)
    if isinstance(data_obj.mo_coeff, (list, tuple)):
        # UHF/UKS: Select coefficients for the specific spin channel
        target_coeffs = data_obj.mo_coeff[spin_idx][:, orbital_index]
    else:
        # RHF/RKS: Only one set of coefficients exists
        target_coeffs = data_obj.mo_coeff[:, orbital_index]

    # MO = Linear combination of AO (Dot product)
    mo_values = np.dot(ao_values, target_coeffs)

    # prepare data for PyVista
    grid_values = mo_values.reshape(nx, ny, nz)

    # create structured Grid
    X, Y, Z = np.meshgrid(x, y, z)
    grid = pv.StructuredGrid(X, Y, Z)

    # assign grid values
    grid_values = mo_values.reshape((nx, ny, nz), order='F')
    grid_values = np.transpose(grid_values, (0, 2, 1)) 
    grid.point_data["values"] = grid_values.flatten()

    # create two Mesh object and assign names
    iso_pos = grid.contour([iso_level], scalars="values")
    iso_neg = grid.contour([-iso_level], scalars="values")

    visual_objects.append((iso_pos))
    visual_objects.append((iso_neg)) 
       
    return visual_objects

def mo_batch(new_data,orb_index,iso_level, nx=50, ny = 50, nz = 50, padding = 3, spin_idx=0):
    visual_objects_raw = draw_orb_molden(new_data, orbital_index=orb_index,spin_idx=spin_idx,
                    iso_level=iso_level, nx = nx, ny = ny, nz = nz, padding = padding)  
    orb_mesh = []
    if hasattr(visual_objects_raw, "n_blocks"):
    # convert MultiBlock object into mesh list
        visual_objects = [visual_objects_raw[i] for i in range(visual_objects_raw.n_blocks)]
    else:
    # it is already a list (or iterable object)
        visual_objects = list(visual_objects_raw)
    for i, item in enumerate(visual_objects):
        # in case draw_orb_molden returns Tupel (mesh, args)
        if isinstance(item, tuple):
            mesh, args = item
        else:
        # in case it returns Mesh directly
            mesh = item
            args = {}
        orb_mesh.append(mesh)
    mesh1 = orb_mesh[0]
    mesh2 = orb_mesh[1]  
    return mesh1, mesh2

def get_optimal_cores():
    # 1. macOS (Apple Silicon Check)
    if platform.system() == "Darwin":
        try:
            # Versuche die Performance-Kerne zu finden
            output = subprocess.check_output(['sysctl', '-n', 'hw.perflevel0.physicalcpu'], 
                                           stderr=subprocess.DEVNULL) 
            return max(1, int(output.strip()))
        except Exception:
            # in case sysctl fails
            pass

    # 2. Fallback for Windows, Linux and Intel-Macs
    # physical cores minus 1
    cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
    return max(1, cores - 1)

# stabilize orbital phase during animation!!
def fix_orbital_phase(current_coeffs, previous_coeffs):
# range(current_coeffs.shape[1]) number for orbitals
    overlaps = np.sum(current_coeffs * previous_coeffs, axis=0)
    current_coeffs[:, overlaps < 0] *= -1.0
    return current_coeffs

def get_validated_dm(mol, mo_coeff, mo_occ):
    # 1. UKS (Radical/Tupel) vs. RKS (Array)
    if isinstance(mo_occ, (list, tuple)):
        # mo_coeff[0] = Alpha-Orbitals, mo_occ[0] = Alpha-Occupation
        dm_a = scf.hf.make_rdm1(mo_coeff[0], mo_occ[0])
        dm_b = scf.hf.make_rdm1(mo_coeff[1], mo_occ[1])
        dm = dm_a + dm_b
        target_n = np.sum(mo_occ[0]) + np.sum(mo_occ[1])
    else:
        # RKS/RHF (Simple Array)
        dm = scf.hf.make_rdm1(mo_coeff, mo_occ)
        target_n = np.sum(mo_occ)
    
    # 2. Normalization (for Molden-Import)
    s = mol.intor('int1e_ovlp')
    current_n = np.trace(dm @ s)
    if not np.isclose(current_n, target_n, atol=1e-5):
        dm *= (target_n / current_n)
    return dm

def draw_esp_molden(data_obj, iso_val=0.001, nx = 50, ny = 50, nz = 50, padding = 4):
    mol = data_obj.mol
    to_ang = 0.529177249
    padding_ang = padding

    # 1. use Bohr for calculation
    at_coords_bohr = mol.atom_coords()
    pad_bohr = padding_ang / to_ang # convert padding
    
    mins = np.min(at_coords_bohr, axis=0) - pad_bohr
    maxs = np.max(at_coords_bohr, axis=0) + pad_bohr
    
    # 2. Grid Axes in Bohr
    x_b = np.linspace(mins[0], maxs[0], nx)
    y_b = np.linspace(mins[1], maxs[1], ny)
    z_b = np.linspace(mins[2], maxs[2], nz)
    
    # 3. Calculate Density (PySCF uses Bohr)
    coords_bohr = lib.cartesian_prod([x_b, y_b, z_b])
    ao_values = mol.eval_gto("GTOval", coords_bohr)
    dm = get_validated_dm(mol, data_obj.mo_coeff, data_obj.mo_occ)
    rho = dft.numint.NumInt().eval_rho(mol, ao_values, dm)
    
    # 4. PyVista Grid in Ångström 
    X, Y, Z = np.meshgrid(x_b * to_ang, y_b * to_ang, z_b * to_ang, indexing='ij')
    grid = pv.StructuredGrid(X, Y, Z)
    grid.point_data["rho"] = rho.reshape((nx, ny, nz), order='F').ravel(order='C')
    surf = grid.contour([iso_val], scalars="rho")
    surf = surf.decimate(0.01) 

    points_bohr = surf.points / to_ang
    
    # Electron Block
    # int1e_grids uses multiprocessing via C
    v_mat = mol.intor('int1e_grids', grids=points_bohr) # Form: (n_pts, n_ao, n_ao)
    # Extremely fast contraction via einsum 
    v_elec = -np.einsum('kij,ij->k', v_mat, dm)
    del v_mat # set RAM free
    # Nuclear Block (vectorized)
    at_coords = mol.atom_coords()
    at_charges = mol.atom_charges()
    # fast distance matrix w/o explicit loops
    diff = points_bohr[:, None, :] - at_coords[None, :, :]
    dists = np.sqrt(np.einsum('ijk,ijk->ij', diff, diff))
    v_nuc = np.sum(at_charges / dists, axis=1)
    # ESP
    total_esp = v_elec + v_nuc

    surf.point_data["ESP"] = total_esp

    return surf, "ESP"

def draw_spin(data_obj, iso_val=0.02, nx = 50, ny = 50, nz = 50, padding = 4):
    mol, mo_coeff, mo_occ = data_obj.mol, data_obj.mo_coeff, data_obj.mo_occ
    
    # Check if we have an unrestricted system (tuple/list for alpha and beta)
    if isinstance(mo_coeff, (list, tuple)):
        # UHF/UKS: Generate two density matrices and sum them up
        dm_alpha = scf.hf.make_rdm1(mo_coeff[0], mo_occ[0])
        dm_beta = scf.hf.make_rdm1(mo_coeff[1], mo_occ[1])
         # Re-normalization (use actual occupation)
        s = mol.intor('int1e_ovlp')
        target_a = np.sum(mo_occ[0])
        target_b = np.sum(mo_occ[1])
        
        # Apply correction for each spin 
        dm_alpha *= (target_a / np.trace(dm_alpha @ s))
        dm_beta *= (target_b / np.trace(dm_beta @ s))
    else:
        # RHF/RKS: Generate standard density matrix
        dm_alpha=scf.hf.make_rdm1(mo_coeff, mo_occ)
        dm_beta=dm_alpha

    to_ang = 0.529177249
    padding_ang = padding

    # 1. Calculation in Bohr 
    at_coords_bohr = mol.atom_coords()
    pad_bohr = padding_ang / to_ang # convert padding
    
    mins = np.min(at_coords_bohr, axis=0) - pad_bohr
    maxs = np.max(at_coords_bohr, axis=0) + pad_bohr
    
    # 2. Grid Axes in Bohr
    x_b = np.linspace(mins[0], maxs[0], nx)
    y_b = np.linspace(mins[1], maxs[1], ny)
    z_b = np.linspace(mins[2], maxs[2], nz)
    
    # 3. Calculate Density (PySCF uses Bohr)
    coords_bohr = lib.cartesian_prod([x_b, y_b, z_b])
    ao_values = mol.eval_gto("GTOval", coords_bohr)
    # --- Density ---
    rho_alpha = dft.numint.NumInt().eval_rho(mol, ao_values, dm_alpha)
    rho_beta = dft.numint.NumInt().eval_rho(mol, ao_values, dm_beta)
    rho_spin = rho_alpha -  rho_beta

    # 4. PyVista Grid in Ångström 
    X, Y, Z = np.meshgrid(x_b * to_ang, y_b * to_ang, z_b * to_ang, indexing='ij')
    grid = pv.StructuredGrid(X, Y, Z)
    grid.point_data["spin"] = rho_spin.reshape((nx, ny, nz), order='F').ravel(order='C')

    # create two Mesh object and assign names
    iso_pos = grid.contour([iso_val], scalars="spin")
    iso_neg = grid.contour([-iso_val], scalars="spin")

    visual_objects = [] # list object for mesh and color
    visual_objects.append((iso_pos))
    visual_objects.append((iso_neg)) 
       
    return visual_objects
    
def draw_spin_mapped(data_obj, iso_val=0.002, nx = 50, ny = 50, nz = 50, padding = 4):
    mol, mo_coeff, mo_occ = data_obj.mol, data_obj.mo_coeff, data_obj.mo_occ
    
    # Check if we have an unrestricted system (tuple/list for alpha and beta)
    if isinstance(mo_coeff, (list, tuple)):
        # UHF/UKS: Generate two density matrices and sum them up
        dm_alpha = scf.hf.make_rdm1(mo_coeff[0], mo_occ[0])
        dm_beta = scf.hf.make_rdm1(mo_coeff[1], mo_occ[1])
         # 2. Renormierung auf die physikalischen Besetzungszahlen
        s = mol.intor('int1e_ovlp')
        target_a = np.sum(mo_occ[0])
        target_b = np.sum(mo_occ[1])
        
        # Korrekturfaktor für jeden Kanal einzeln (wegen des Basis-Mismatchs)
        dm_alpha *= (target_a / np.trace(dm_alpha @ s))
        dm_beta *= (target_b / np.trace(dm_beta @ s))
    else:
        # RHF/RKS: Generate standard density matrix
        dm_alpha=scf.hf.make_rdm1(mo_coeff, mo_occ)
        dm_beta=dm_alpha

    to_ang = 0.529177249
    padding_ang = padding
    
    # 1. Calculation in Bohr 
    at_coords_bohr = mol.atom_coords()
    pad_bohr = padding_ang / to_ang # convert padding 
    
    mins = np.min(at_coords_bohr, axis=0) - pad_bohr
    maxs = np.max(at_coords_bohr, axis=0) + pad_bohr
    
    # 2. Grid Axes in Bohr
    x_b = np.linspace(mins[0], maxs[0], nx)
    y_b = np.linspace(mins[1], maxs[1], ny)
    z_b = np.linspace(mins[2], maxs[2], nz)
    
    # 3. Calculate Density (PySCF uses Bohr)
    coords_bohr = lib.cartesian_prod([x_b, y_b, z_b])
    ao_values = mol.eval_gto("GTOval", coords_bohr)
    # --- Density ---
    rho_alpha = dft.numint.NumInt().eval_rho(mol, ao_values, dm_alpha)
    rho_beta = dft.numint.NumInt().eval_rho(mol, ao_values, dm_beta)

    rho_spin = rho_alpha - rho_beta
    rho_total = rho_alpha + rho_beta # physical unit: e/Bohr^3

    # 4. PyVista Grid in Ångström 
    X, Y, Z = np.meshgrid(x_b * to_ang, y_b * to_ang, z_b * to_ang, indexing='ij')
    grid = pv.StructuredGrid(X, Y, Z)
    grid.point_data["spin_abs"] = rho_spin.reshape((nx, ny, nz), order='F').ravel(order='C')
    grid.point_data["total_rho"] = rho_total.reshape((nx, ny, nz), order='F').ravel(order='C')

    # create geometry based on physical density
    surf = grid.contour([iso_val], scalars="total_rho")
    surf = surf.interpolate(grid)

    # Polarization: in percentage of alpha-spin
    surf.point_data["mapped_data"] = surf.point_data["spin_abs"] / iso_val

    return surf, "mapped_data"

################# POV-Ray Export ###################################
def export_pov_mo(mesh1, mesh2, filename="test.inc", object_name="name", idx=0):
    #MeshObject
    # prepare 1. Mesh
    if mesh1.n_points > 0 and mesh1.n_cells > 0:
        mesh1 = mesh1.extract_surface(algorithm='dataset_surface').triangulate().compute_normals(cell_normals=False, point_normals=True)
        normals = mesh1.point_data["Normals"]
        verts = mesh1.points
        faces = mesh1.faces.reshape(-1, 4)[:, 1:]

        with open(filename, 'a') as f:
            f.write(f"#declare MeshObject_1_{idx} = mesh2 {{\n")
        
            # --- A. Vertex-Vector (coordinates) 
            f.write("  vertex_vectors {\n")
            f.write(f"    {len(verts)},\n")
            for v in verts:
                f.write(f"    <{v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f}>,\n")
            f.write("  }\n")

            # --- B. Normal Vectors ---
            f.write("  normal_vectors {\n")
            f.write(f"    {len(normals)},\n")
            for n in normals:
                f.write(f"    <{n[0]:.6f}, {n[1]:.6f}, {n[2]:.6f}>,\n")
            f.write("  }\n")

            # --- C. Point Indices (Geometry) ---
            f.write("  face_indices {\n")
            f.write(f"    {len(faces)},\n")
            for face in faces:
                f.write(f"    <{face[0]}, {face[1]}, {face[2]}>,\n")
            f.write("  }\n")

            # --- D. Normal Indices (Smooting) ---
            f.write("  normal_indices {\n")
            f.write(f"    {len(faces)},\n")
            for face in faces:
                f.write(f"    <{face[0]}, {face[1]}, {face[2]}>,\n")
            f.write("  }\n")
        
            f.write("}\n")
    else:
    # in case of an empty mesh 
        with open(filename, 'a') as f:
            f.write(f"#declare MeshObject_1_{idx} = union {{ }}\n")

    # prepare 2. Mesh 
    if mesh2.n_points > 0 and mesh2.n_cells > 0:
        mesh2 = mesh2.extract_surface(algorithm='dataset_surface').triangulate().compute_normals(cell_normals=False, point_normals=True)
        normals = mesh2.point_data["Normals"]
        verts = mesh2.points
        faces = mesh2.faces.reshape(-1, 4)[:, 1:]

        with open(filename, 'a') as f:
            f.write(f"#declare MeshObject_2_{idx} = mesh2 {{\n")
            
            # --- A. Vertex-Vectors (coordinates)
            f.write("  vertex_vectors {\n")
            f.write(f"    {len(verts)},\n")
            for v in verts:
                f.write(f"    <{v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f}>,\n")
            f.write("  }\n")

            # --- B. Normal-Vectors ---
            f.write("  normal_vectors {\n")
            f.write(f"    {len(normals)},\n")
            for n in normals:
                f.write(f"    <{n[0]:.6f}, {n[1]:.6f}, {n[2]:.6f}>,\n")
            f.write("  }\n")

            # --- C. Point-Indices (Geometry) ---
            f.write("  face_indices {\n")
            f.write(f"    {len(faces)},\n")
            for face in faces:
                f.write(f"    <{face[0]}, {face[1]}, {face[2]}>,\n")
            f.write("  }\n")

            # --- D. Normal-Indices (Smoothing) ---
            f.write("  normal_indices {\n")
            f.write(f"    {len(faces)},\n")
            for face in faces:
                f.write(f"    <{face[0]}, {face[1]}, {face[2]}>,\n")
            f.write("  }\n")
            
            f.write("}\n")
    else:
        # in case of an empty mesh
        with open(filename, 'a') as f:
            f.write(f"#declare MeshObject_2_{idx} = union {{ }}\n")
 
# combine all objects
    with open(filename, 'a') as f:
        f.write(f"""\
// Combine Atoms, Bonds and Mesh 
#declare {object_name}[{idx}] = union {{
    object {{AtomsGroup_{idx}}}
    object {{BondsGroup_{idx}}}
    object {{ MeshObject_1_{idx} texture {{ pigment {{ color color_pos filter trans }} finish {{OrbFinish}} }} }}
    object {{ MeshObject_2_{idx} texture {{ pigment {{ color color_neg filter trans }} finish {{OrbFinish}} }} }}
}}
// End of Section {idx}
        """)

def export_pov_esp(mesh,filename="test.inc", object_name="name", 
                   cmap_name="bwr", clim=[-0.02,0.02], idx=0):
    # 1. prepare Mesh Object
    mesh=mesh.mapper.dataset
    mesh = mesh.triangulate()
    verts = mesh.points
    faces = mesh.faces.reshape(-1, 4)[:, 1:]
    
    # 2. calculate colors (mapping of ESP values)
    # esp values in mesh.active_scalars 
    scalars = mesh.active_scalars
    norm = plt.Normalize(vmin=clim[0], vmax=clim[1])
    colormap = colormaps.get_cmap(cmap_name)
    
    # RGB colors for each Vertex (0.0 to 1.0 for POV-Ray)
    colors = colormap(norm(scalars))[:, :3] 

    # Calculate Normals in PyVista
    mesh = mesh.compute_normals(cell_normals=False, point_normals=True)
    normals = mesh.point_data["Normals"]

    with open(filename, 'a') as f:
        f.write(f"#declare MeshObject_{idx} = mesh2 {{\n")
        
        # Vertices
        f.write("  vertex_vectors {\n")
        f.write(f"    {len(verts)},\n")
        for v in verts:
            f.write(f"    <{v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f}>,\n")
        f.write("  }\n")

        # Normal Vectors
        f.write("  normal_vectors {\n")
        f.write(f"    {len(normals)},\n")
        for n in normals:
            f.write(f"    <{n[0]:.6f}, {n[1]:.6f}, {n[2]:.6f}>,\n")
        f.write("  }\n")
        
        # Textures (Color per Vertex)
        f.write("  texture_list {\n")
        f.write(f"    {len(verts)},\n")
        for c in colors:
            # define Texture for each point
            f.write(f"    texture {{ pigment {{ color rgb <{c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f}> filter trans }}  }}\n")
        f.write("  }\n")
        
        # Faces with Texture Interpolation
        f.write("  face_indices {\n")
        f.write(f"    {len(faces)},\n")
        for face in faces:
            # the three numbers following the Mesh Index <v1, v2, v3> 
            # are Textures for each Vertex, 
            # POV-Ray interpolates the color between Vertices
            f.write(f"    <{face[0]}, {face[1]}, {face[2]}>, {face[0]}, {face[1]}, {face[2]},\n")
        f.write("  }\n")
        
        f.write("}\n")

# combine all objects
    with open(filename, 'a') as f:
        f.write(f"""\
// Combine Atoms, Bonds and Mesh 
#declare {object_name}[{idx}] = union {{
    object {{AtomsGroup_{idx}}}
    object {{BondsGroup_{idx}}}
    object {{MeshObject_{idx}}}
}}
// End of Section {idx}
        """)

def export_pov_mol(points, atom_types,cpk_colors=None,filename="test.inc", idx=0):
    #AtomsGroup
    with open(filename, 'a') as f:
        f.write(f"//Begin of Section {idx}  \n")
        f.write(f"#declare AtomsGroup_{idx} = union {{\n")
        for i, pos in enumerate(points):
            val = int(atom_types[i])
            # get colors from Color Dictionary
            color_name = cpk_colors.get(val, "magenta")
            # Conversion of Names to RGB for POV-Ray 
            rgb = {"white": "<1,1,1>", "gray": "<.3,.3,.3>", "blue": "<0,0,1>", 
                "red": "<1,0,0>", "orange": "<1, 0.55, 0>", "yellow": "<1,1,0>", 
                "brown": "<1,0.65,0>", "darkred": "<0.5,0,0>", "green": "<0, 0.82,0>"}.get(color_name, "<1,0,1>")

            atomic_num = int(atom_types[i])
            match atomic_num:
                case 1: # Hydrogen
                    rad_var = "atom_rad_h"
                case _ if 3 <= atomic_num <= 10: # 2nd period (He-Ne)
                    rad_var = "atom_rad_2"
                case _ if 11 <= atomic_num <= 18: # 3rd period (Na-Ar)
                    rad_var = "atom_rad_3"
                case _: # every other element
                    rad_var = "atom_rad_def"

            f.write(f"  sphere {{ <{pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}>, {rad_var}\n")
            f.write(f"    pigment {{ color rgb {rgb} filter trans_atom }}\n")
            f.write("    finish { AtomFinish }\n")
            f.write("  }\n")
        f.write("}\n")
    #BondsGroup
    with open(filename, 'a') as f:
        f.write(f"#declare BondsGroup_{idx} = union {{\n")
        for i in range(len(points)):
            type_i = int(atom_types[i])
            rad_i = cov_radii.get(type_i, default_radius)
            for j in range(i + 1, len(points)):
                type_j = int(atom_types[j])
                rad_j = cov_radii.get(type_j, default_radius) 
                bd_threshold = rad_i + rad_j + 0.6
                dist = np.linalg.norm(points[i] - points[j])
                # Threshold for Bonds in Angstrom
                if 0.6 < dist < bd_threshold:
                    p1 = points[i]
                    p2 = points[j]
                    f.write(f"  cylinder {{ <{p1[0]:.4f}, {p1[1]:.4f}, {p1[2]:.4f}>, "
                            f"<{p2[0]:.4f}, {p2[1]:.4f}, {p2[2]:.4f}>, bond_rad\n")
                    f.write("    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }\n")
                    f.write("    finish { BdFinish }\n")
                    f.write("  }\n")
        f.write("}\n")

def header_mo(rgb_pos,rgb_neg,o_file, obj_name, length,trans):
    with open(o_file, 'w') as f:
        f.write(f"""\
// ---- created with BatchMol {ver_no} by Dr. Tobias Schulz
// ---- MO Cube Object Object: #include "{o_file}" in povray
// ---- use "object{{{obj_name}[i]}}" in code
//declare molecule object array
#declare {obj_name} = array[{length+1}];
// ---- declare Variables
//
// ---- Orbital Mesh Section
// transparency and color
#declare trans = {trans}; 
#declare color_pos = rgb <{rgb_pos[0]:.3f}, {rgb_pos[1]:.3f}, {rgb_pos[2]:.3f}>;
#declare color_neg = rgb <{rgb_neg[0]:.3f}, {rgb_neg[1]:.3f}, {rgb_neg[2]:.3f}>;
// 1. defined presets for orbital finishes
#declare Fin_Std      = finish {{ phong 0.3 ambient 0.2 diffuse 0.6 }}
#declare Fin_Glassy   = finish {{ phong 0.9 specular 0.8 reflection 0.1 roughness 0.001 }}
#declare Fin_Metallic = finish {{ phong 0.5 metallic 0.7 brilliance 2.0 diffuse 0.3 }}
#declare Fin_Matte    = finish {{ phong 0.0 ambient 0.1 diffuse 0.8 }}
// 2. select active orbital finish
#declare OrbFinish = Fin_Glassy; 
        """)

def header_esp(o_file,obj_name, length,trans):
    with open(o_file, 'w') as f:
        f.write(f"""\
// ---- created with BatchMol {ver_no} by Dr. Tobias Schulz
// ----  MO Cube Object Object: #include "{o_file}" in povray
// ---- use "object{{{obj_name}[i]}}" in code
// Scalar_Bar is also available as "Colorbar"
//declare molecule object array
#declare {obj_name} = array[{length+1}];
// ---- declare Variables
//
//ESP Mesh Section Transparency
#declare trans = {trans}; 
        """)

def header_spin(rgb_pos,rgb_neg,filename, object_name, length,trans):
    with open(filename, 'w') as f:
        f.write(f"""\
// ---- created with BatchMol {ver_no} by Dr. Tobias Schulz
// Spin Density Object Object: #include "{object_name}.inc" into povray
//use "object{{{object_name}}}" in code
//declare molecule object array
#declare {object_name} = array[{length+1}];
// ---- declare Variables
//
//declare Variables
#declare trans = {trans}; //Mesh Transparency
#declare color_pos = rgb <{rgb_pos[0]:.3f}, {rgb_pos[1]:.3f}, {rgb_pos[2]:.3f}>;
#declare color_neg = rgb <{rgb_neg[0]:.3f}, {rgb_neg[1]:.3f}, {rgb_neg[2]:.3f}>;
// Pre-defined Finishes
// 1. Die Presets definieren (mit umschließenden finish-Klammern)
#declare Fin_Std      = finish {{ phong 0.3 ambient 0.2 diffuse 0.6 }}
#declare Fin_Glassy   = finish {{ phong 0.9 specular 0.8 reflection 0.1 roughness 0.001 }}
#declare Fin_Metallic = finish {{ phong 0.5 metallic 0.7 brilliance 2.0 diffuse 0.3 }}
#declare Fin_Matte    = finish {{ phong 0.0 ambient 0.1 diffuse 0.8 }}
// 2. Die Auswahl treffen (das macht dein Python-Export)
#declare OrbFinish = Fin_Glassy; 
//
        """) 

def header_spin_mapped(filename, object_name, length, trans=0.66):
    with open(filename, 'w') as f:
        f.write(f"""\
// ---- created with BatchMol {ver_no} by Dr. Tobias Schulz
// Mapped Spin Density Object: #include "{object_name}.inc" into povray
// use "object{{{object_name}}}" in code
// Scalar_Bar is also available as "Colorbar"
//declare molecule object array
#declare {object_name} = array[{length+1}];
//---- declare Variables
//
#declare trans = {trans}; //Mesh Transparency
//
        """)

def header_mol(o_file):
    with open(o_file, 'a') as f:
        f.write(f"""\
//
// ---- Atom and Bond Section
//transparency
#declare trans_bd = 0;
#declare trans_atom = 0;
//atom radius
#declare atom_rad_h = 0.24;
#declare atom_rad_2 = 0.35;
#declare atom_rad_3 = 0.42;
#declare atom_rad_def = 0.5;
#declare bond_rad = 0.08;
// predefined finishes:
#declare Fin_Glassy   = finish {{ phong 0.9 specular 0.8 reflection 0.1 roughness 0.001 }}
#declare Fin_Metallic = finish {{ phong 0.5 metallic 0.7 brilliance 2.0 diffuse 0.3 }}
#declare Fin_Matte    = finish {{ phong 0.0 ambient 0.1 diffuse 0.8 }}
// Define Bond Finishes
#declare Fin_Bd_Std = finish {{ phong 0.2 ambient 0.2 }}
// Select bond Finish
#declare BdFinish = Fin_Bd_Std;
// Definde Atom finishes
#declare Fin_Atom_Std = finish {{ phong 0.6 specular 0.4 ambient 0.2 }}
// select Atom Finish
#declare AtomFinish = Fin_Atom_Std;
        """)

def export_pov_colorbar(filename, cmap_name, clim, type, height=2.0, radius=0.08):
    """
    Creates Cylinder/Legend for ColorMapping v_min/v_max.
    """
    v_min, v_max = clim
    cmap = plt.get_cmap(cmap_name)
    num_samples = 10 
    
    # positional argument for text label
    text_offset_x = radius + 0.2
    font_name = "arial.ttf" # font must exist in POV-Ray path

    with open(filename, 'a') as f:
        f.write("\n// --- Colorbar with Labels ---\n")
        # 1. Cylinder
        f.write(f"#declare Bar = cylinder {{ <0, 0, 0>, <0, {height}, 0>, {radius}\n")
        f.write("  pigment { gradient y color_map {\n")
        for i in range(num_samples):
            frac = i / (num_samples - 1)
            rgba = cmap(frac)
            f.write(f"    [{frac:.3f} color rgb <{rgba[0]:.4f}, {rgba[1]:.4f}, {rgba[2]:.4f}>]\n")
        f.write(f"  }} scale <1, {height}, 1> }}\n")
        f.write("  finish { ambient 0.7 diffuse 0.3 }\n}\n")

        # 2. Text Label
        # text { ttf "font.ttf" "string" thickness, offset }
        label_min = f'text {{ ttf "{font_name}" "{v_min:.2f}" 0.05, 0 pigment {{ rgb 1 }} finish {{ ambient 1 diffuse 0 }} scale 0.35 }}'
        label_max = f'text {{ ttf "{font_name}" "{v_max:.2f}" 0.05, 0 pigment {{ rgb 1 }} finish {{ ambient 1 diffuse 0 }} scale 0.35 }}'
        if type == "esp": title = "ESP(Hartree)"
        else: title = "Spindensity (au)"
        label_legend = f'text {{ ttf "{font_name}" "{title}" 0.05, 0 pigment {{ rgb 1 }} finish {{ ambient 1 diffuse 0 }} scale 0.35 }}'

        # 3. Combine Cylinder and Text Label
        f.write("#declare Colorbar = union {\n")
        f.write("  object { Bar }\n")
        f.write(f"  object {{ {label_min} translate <{text_offset_x}, 0, 0> }}\n")
        f.write(f"  object {{ {label_max} translate <{text_offset_x}, {height - 0.2}, 0> }}\n")
        f.write(f"  object {{ {label_legend} translate <0, {height + 0.4}, 0> }}\n")
        f.write("}\n")

#### BLENDER ColorBar/Scalebar ######
def create_3d_colorbar_group(v_min, v_max, mode="esp", cmap_name="rainbow", height=1.5, width=0.15):
    """
    Creates a compact 3D colorbar. 
    Height reduced from 3.0 to 1.5, width from 0.3 to 0.15.
    """
    visuals = [] 
    
    # 1. Rectangular Bar
    # Position shifted closer to the center (x=3.5 instead of 5)
    cb_poly = pv.Plane(center=(3.5, 0.0, 0.0), direction=(1, 0, 0), 
                       i_size=width, j_size=height)
    
    y_min, y_max = np.min(cb_poly.points[:, 1]), np.max(cb_poly.points[:, 1])
    norm_values = np.interp(cb_poly.points[:, 1], (y_min, y_max), (v_min, v_max))
    cb_poly.point_data["cb_scalars"] = norm_values
    
    # Material name 'emit_scale' helps us identify it in Blender for emission setup
    visuals.append((cb_poly, "scale_bar_emit", {
        "scalars": "cb_scalars", 
        "cmap": cmap_name, 
        "clim": [v_min, v_max]
    }))

    # 2. Labels (Smaller text)
    label_args = {"color": "white", "smooth_shading": True}
    text_scale = 0.12  # Reduced from 0.25
    text_depth = 0.02
    
    # v_min Label (bottom)
    txt_min = pv.Text3D(f"{v_min:.2f}", depth=text_depth)
    txt_min.scale(text_scale)
    # Positions adjusted for smaller height
    txt_min.translate([3.7, -height/2, 0])
    visuals.append((txt_min, "lbl_min_emit", label_args))
    
    # v_max Label (top)
    txt_max = pv.Text3D(f"{v_max:.2f}", depth=text_depth)
    txt_max.scale(text_scale)
    txt_max.translate([3.7, height/2 - 0.1, 0])
    visuals.append((txt_max, "lbl_max_emit", label_args))
    
    # Title (above)
    title_str = "ESP [a.u.]" if mode == "esp" else "Spin [a.u.]"
    txt_title = pv.Text3D(title_str, depth=text_depth)
    txt_title.scale(text_scale * 1.2)
    txt_title.translate([3.4, height/2 + 0.2, 0])
    visuals.append((txt_title, "lbl_title_emit", label_args))
    
    return visuals

### BLENDER Script #####
def generate_blender_script(path):
    """
    Generates a companion Blender Python script for Multi-File GLB export.
    Frozen State: Correct sorting under Trajectory_Control, Slot-Logic, 
    and robust Constant Interpolation.
    """
    script_path = os.path.splitext(path)[0] + "_setup.py"
    current_path = os.getcwd().replace("\\", "/") 
    
    blender_script = f"""# created with BatchMol {ver_no} (Multi-File Orbitals Frozen)
# ==============================================================================
# USER GUIDE: 
# 1. RUN THIS SCRIPT
# 2. CLEANUP: Search 'Renderer Node' in Outliner -> Select all (A) 
#    -> Right Click -> 'Delete Hierarchy'. This removes Cameras but keeps Meshes.
# ==============================================================================

import bpy, os, re

# --- Settings ---
path_to_glb = "{current_path}"
protected = ["Camera", "Plane", "Cylinder", "Sun", "World", "TRAJECTORY_CONTROL", "MASTER", "DUMMY"]

# 1. Setup Controller
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. File Discovery
files = sorted([f for f in os.listdir(path_to_glb) if f.endswith(".glb")])
frame_files = [f for f in files if "scalebar" not in f]

# 3. Import Loop
for i, filename in enumerate(frame_files):
    filepath = os.path.join(path_to_glb, filename)
    bpy.ops.import_scene.gltf(filepath=filepath)
    
    # Track imported objects
    new_objs = [o for o in bpy.context.selected_objects]
    # Sort them by their internal Blender name (mesh_0, mesh_1...)
    new_objs.sort(key=lambda o: [int(c) if c.isdigit() else c.lower() for c in re.split('(\\\\d+)', o.name)])
    
    current_frame = i + 1 

    # We only care about Meshes for materials and animation
    meshes_in_file = [o for o in new_objs if o.type == 'MESH']
    
    for slot_idx, obj in enumerate(meshes_in_file):
        # A) Material Assignment (Slot Logic)
        m_names = ["MASTER_Molecule", "MASTER_Orb_Pos", "MASTER_Orb_Neg", "MASTER_Surface"]
        if slot_idx < len(m_names):
            mat = bpy.data.materials.get(m_names[slot_idx])
            if mat:
                obj.data.materials.clear()
                obj.data.materials.append(mat)
                obj.color = mat.diffuse_color

        # B) Parenting & Animation
        obj.parent = cntrl
        
        # Keyframe Sequence
        obj.scale = (0,0,0)
        obj.keyframe_insert(data_path="scale", frame=current_frame - 1)
        
        # Show only real geometry (Dummies stay hidden)
        if len(obj.data.vertices) > 1:
            obj.scale = (1,1,1)
        obj.keyframe_insert(data_path="scale", frame=current_frame)
        
        obj.scale = (0,0,0)
        obj.keyframe_insert(data_path="scale", frame=current_frame + 1)
        
        # C) Robust Constant Interpolation
        if obj.animation_data and obj.animation_data.action:
            action = obj.animation_data.action
            if hasattr(action, "fcurves"):
                for fc in action.fcurves:
                    for kp in fc.keyframe_points:
                        kp.interpolation = 'CONSTANT'

    # 4. Cleanup Hierarchy: Remove empty containers (keeps Meshes due to parenting)
    for o in new_objs:
        if o.name in bpy.data.objects:
            if o.type == 'EMPTY' or not o.data:
                bpy.data.objects.remove(o, do_unlink=True)

# Finalize
bpy.context.scene.frame_start = 1
bpy.context.scene.frame_end = len(frame_files)
bpy.context.scene.frame_set(1)

print(f"Multi-File Sync finished: {{len(frame_files)}} frames ready.")
"""
    with open(script_path, "w") as f:
        f.write(blender_script)
    print(f"Frozen Multi-File Script written to: {script_path}")

def generate_blender_script_one(path):
    script_path = os.path.splitext(path)[0] + "_setup.py"
    blender_script = f"""# created with BatchMol {ver_no} (One-File Orbitals Frozen)
import bpy, re, os

# --- Helper: Natural Sort (mesh_2 comes before mesh_10) ---
def natural_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split('(\\\\d+)', text)]

# 1. Controller Setup
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. Discovery and FORCE SORT
container = bpy.data.objects.get("Renderer Node")
if not container:
    all_objs = [o for o in bpy.data.objects if o.type == 'MESH' and "MASTER" not in o.name]
else:
    all_objs = [child for child in container.children if child.type == 'MESH']

# HIER PASSIERT DIE MAGIE: Wir zwingen Blender in die richtige Reihenfolge
all_objs.sort(key=lambda o: natural_key(o.name))

# 3. Processing Slots
slots_per_frame = 4
for i, obj in enumerate(all_objs):
    frame_idx = (i // slots_per_frame) + 1
    slot_idx = i % slots_per_frame 
    
    obj.parent = cntrl
    
    # Material Assignment
    m_names = ["MASTER_Molecule", "MASTER_Orb_Pos", "MASTER_Orb_Neg", "MASTER_Surface"]
    mat = bpy.data.materials.get(m_names[slot_idx])
    if mat:
        obj.data.materials.clear()
        obj.data.materials.append(mat)
        obj.color = mat.diffuse_color

    # Animation
    is_dummy = len(obj.data.polygons) == 0
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx - 1)
    obj.scale = (1, 1, 1) if not is_dummy else (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx)
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx + 1)

    # --- Constant Interpolation Fix ---
    if obj.animation_data and obj.animation_data.action:
        act = obj.animation_data.action
        if hasattr(act, "fcurves"):
            for fc in act.fcurves:
                for kp in fc.keyframe_points:
                    kp.interpolation = 'CONSTANT'

# 4. Viewport Fix
for area in bpy.context.screen.areas:
    if area.type == 'VIEW_3D':
        for space in area.spaces:
            if space.type == 'VIEW_3D': space.shading.color_type = 'OBJECT'

bpy.context.scene.frame_end = len(all_objs) // slots_per_frame
bpy.context.scene.frame_set(1)
print(f"Sorted {{len(all_objs)}} objects. Timeline ready.")
"""
    with open(script_path, "w") as f:
        f.write(blender_script)


cpk_colors = defaultdict(lambda: "magenta")
cpk_colors.update({
            1: "white",  #H
            5: "pink",  #B
            6: "gray",   #C
            7: "blue",   #N
            8: "red",    #O
            9: "orange",  #F
            14: "darkgrey", #Si
            12: "darkgreen", #Mg
            15: "brown",  #P
            16: "yellow", #S
            17: "green",  #Cl
            26: "darkorange", #Fe
            24: "darkcyan",    # Cr (Chrom)
            27: "royalblue",   # Co (Cobalt)
            28: "silver",      # Ni (Nickel)
            29: "chocolate",   # Cu (Kupfer)
            40: "cadetblue",   # Zr (Zirconium)
            44: "teal",        # Ru (Ruthenium)
            45: "deeppink",    # Rh (Rhodium) 
            78: "lightgrey",    # Pt (Platin)
            35: "darkred",   #Br
            53: "darkviolet" # I
        })
        # Bond Parameters
cov_radii = {
        1: 0.31,   # H
        5: 0.82,   # B
        6: 0.76,   # C
        7: 0.71,   # N
        8: 0.66,   # O
        9: 0.57,   # F
        14: 1.11,  # Si
        15: 1.06,  # P
        16: 1.05,  # S
        17: 1.02,  # Cl
        35: 1.20,  # Br
        53: 1.39,  # I
        24: 1.39,  # Cr
        27: 1.26,  # Co
        28: 1.21,  # Ni
        29: 1.32,  # Cu
        40: 1.48,  # Zr
        44: 1.26,  # Ru
        45: 1.35,  # Rh
        78: 1.28   # Pt
        }
        # Standard Radius for unknown elements
default_radius = 1.0

# current Version !!!!
ver_no = "4.3"
# new:
# - corrected ESP and Spin Density
# - adaptive parallelized mode for ESP calculation
# - stabilized color dictionary
# - blender: adaptive atom radius and bond threshold
# - povray: adaptive bond threshold
# - improved parallelization, requires parallel PySCF
# - improved iso density calculation
# - Scalebar for Blender export

@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument('files', nargs=-1) # 'e.g. scan.*.molden'
@click.option('--o_file', '-o', default='output',
                help='output file, default output.inc/.glb') # "scan.inc"
@click.option('--obj_name', '-n', default='molecule',
              help='object name in output.inc, default "molecule" ') #scan
@click.option('--iso_level', '-l', type=float,
              help='iso level, default 0.02 (orbitals), 0.002 (iso-surface)')
@click.option('--n_pts', nargs = 3, type=int, default=(50,50,50),
              help='number of grid points, default 50 50 50')
@click.option('--padding', type=float, default=4,
              help='padding for grid calculation, default 4')
@click.option('--type', '-p', default='esp',
              help='execution mode "mo", "esp", "spin", "spin-m", default "esp" ')
@click.option('--cmap','-c', default='rainbow',
              help='colormap for esp plot, default "rainbow"')
@click.option('--orb_index','-i', type=int, default=0,
              help='orbital index for mo plot, default "0"')
@click.option('--spin','-s', type=int, default=0,
              help='alpha or beta spin default "0" = alpha')
@click.option('--v_max', '-v', type=float, default=None,
              help='limit for scalebar, v_min=-v_max, default +/-0.15 ESP, +/-1 Spin')
@click.option('--o_mode', '-M', default="pov",
              help='output mode "pov": POV-Ray *.inc (default), "bld": Blender Multi *.glb, "bld-one": Blender OneFile *.glb')

def main(files,o_file,obj_name,iso_level, n_pts, padding, type,cmap,orb_index,spin,v_max,o_mode):
    if not files:
        msg="no input found, normal termination"
        click.echo(msg)
        return
    trans=0.6
    i_files=files
    length=len(i_files)

    if v_max is None:
        if type == "spin-m":
            v_max = 1.0
            print(f"Notice: No v_max provided. Using {v_max} for Spin-Mapping.")
        elif type == "esp": 
            v_max = 0.15
            print(f"Notice: No v_max provided. Using {v_max} for ESP-Mapping.")
    v_min = -v_max if v_max is not None else None
    
    if iso_level is None:
        if type in ["spin-m", "esp"]:
            iso_level = 0.002
            print(f"Notice: No iso_level provided. Using {iso_level} for Iso-Surface.")
        elif type in ["mo", "spin"]: 
            iso_level = 0.02
            print(f"Notice: No iso_level provided. Using {iso_level} for Orbitals.")

    #grid variables
    nx, ny, nz = n_pts
    padding = padding
            
    #orbital colors
    color_pos = "#F60505F3"
    color_neg = "#5207F3F3"
    rgb_pos = pv.Color(color_pos).float_rgb
    rgb_neg = pv.Color(color_neg).float_rgb
    if o_mode == "pov":
        o_file += ".inc"
        #prepare header section
        match type:    # export general information and variable declaration
            case "mo":
                header_mo(rgb_pos, rgb_neg,o_file,obj_name,length,trans)
            case "esp":
                header_esp(o_file,obj_name, length,trans)
                export_pov_colorbar(o_file, cmap, [v_min,v_max], type)
            case "spin":
                header_spin(rgb_pos, rgb_neg,o_file,obj_name,length,trans)
            case "spin-m":
                header_spin_mapped(o_file,obj_name, length,trans)
                export_pov_colorbar(o_file, cmap, [v_min,v_max], type)
        header_mol(o_file) # export atom and bond section
        plotter = pv.Plotter(off_screen=True)
        if i_files:
            p_bar = tqdm(i_files)
            old_mo_coeffs = None
            for i, name in enumerate(p_bar, start=1):
                p_bar.set_description(f"processing {name}")
                ext = name.split('.')[-1]
                if ext == 'molden':
                    new_data=MoleculeData.from_molden(name)
                elif ext == 'fchk':
                    new_data=MoleculeData.from_fchk(name)
                else: 
                    print(f"invalid input format: {ext}")
                    sys.exit(1)

                atom_points = new_data.atom_points
                atom_types = new_data.atom_types
            
                match type:
                    case "mo": 
                        if old_mo_coeffs is not None:  
                            if isinstance(new_data.mo_coeff, (list, tuple)):
                            # UHF/UKS: Select coefficients for the specific spin channel
                                fix_orbital_phase(new_data.mo_coeff[spin], old_mo_coeffs[spin])
                            else:
                            # RHF/RKS: Only one set of coefficients exists
                                if new_data.mo_coeff.shape == old_mo_coeffs.shape:
                                    fix_orbital_phase(new_data.mo_coeff, old_mo_coeffs)
                        old_mo_coeffs = np.copy(new_data.mo_coeff) # for comparison of orbital phases
                        
                        mesh1, mesh2 = mo_batch(new_data, orb_index=orb_index,spin_idx=spin, padding= padding, 
                                                iso_level=iso_level, nx = nx, ny = ny, nz = nz)  
                        export_pov_mol(atom_points, atom_types,cpk_colors=cpk_colors,
                                    filename=o_file, idx=i)
                        export_pov_mo(mesh1, mesh2,filename=o_file, 
                            object_name=obj_name, idx=i)
                    case "esp":
                        esp_mesh, active_scalar = draw_esp_molden(new_data,iso_val=iso_level, nx=nx, ny=ny, nz=nz, padding=padding)
                        mesh_args = {
                            "mesh":esp_mesh,
                            "scalars":active_scalar,
                            "cmap":cmap, 
                            "clim":[v_min, v_max], 
                            "opacity":trans,
                            "smooth_shading":True
                            }
                        plotter.clear_actors()
                        ESP_mesh = plotter.add_mesh(**mesh_args)
                        export_pov_mol(atom_points, atom_types,cpk_colors=cpk_colors,
                                    filename=o_file, idx=i)
                        export_pov_esp(ESP_mesh, filename=o_file,object_name=obj_name, 
                                    cmap_name=cmap, clim=[v_min,v_max],idx=i)  
                    case "spin":
                        mesh1, mesh2 = draw_spin(new_data,iso_val=iso_level, nx=nx, ny=ny, nz=nz, padding=padding)
                        export_pov_mol(atom_points, atom_types,cpk_colors=cpk_colors,
                                    filename=o_file, idx=i)
                        export_pov_mo(mesh1, mesh2,filename=o_file, 
                            object_name=obj_name, idx=i)
                    case "spin-m":
                        esp_mesh, active_scalar = draw_spin_mapped(new_data,iso_val=iso_level,nx=nx,ny=ny,nz=nz,padding=padding)
                        mesh_args = {
                            "mesh":esp_mesh,
                            "scalars":active_scalar,
                            "cmap":cmap, 
                            "clim":[v_min, v_max], 
                            "opacity":trans,
                            "smooth_shading":True
                            }
                        plotter.clear_actors()
                        ESP_mesh = plotter.add_mesh(**mesh_args)  
                        export_pov_mol(atom_points, atom_types,cpk_colors=cpk_colors,
                                    filename=o_file, idx=i)
                        export_pov_esp(ESP_mesh, filename=o_file,object_name=obj_name, 
                                    cmap_name=cmap, clim=[v_min,v_max],idx=i) 
            print(f"\nDone! POV-Ray *.inc File written to {o_file}.glb")
    elif o_mode == "bld":
        if type in ["esp", "spin-m"]:
            cb_pl = pv.Plotter()
            cb_group = create_3d_colorbar_group(v_min, v_max, type, cmap)
            for mesh, base_name, kwargs in cb_group:
                cb_pl.add_mesh(mesh, name=base_name, **kwargs)
            cb_pl.export_gltf(f"{o_file}_scalebar.glb")
        if i_files:
            p_bar = tqdm(i_files)
            old_mo_coeffs = None
            for i, name in enumerate(p_bar):  # start i=0
                p_bar.set_description(f"processing {name}")
                new_data=MoleculeData.from_molden(name)
                pl = pv.Plotter() #!!! 
                atom_points = new_data.atom_points
                atom_types = new_data.atom_types

                comb_mesh = draw_mol(atom_points, atom_types, cpk_colors)
                pl.add_mesh(comb_mesh, name=f"mol_{i:03d}", scalars="RGB", rgb=True, smooth_shading=True)

                match type:
                    case "mo":
                        if old_mo_coeffs is not None:  
                            if isinstance(new_data.mo_coeff, (list, tuple)):
                            # UHF/UKS: Select coefficients for the specific spin channel
                                fix_orbital_phase(new_data.mo_coeff[spin], old_mo_coeffs[spin])
                            else:
                            # RHF/RKS: Only one set of coefficients exists
                                if new_data.mo_coeff.shape == old_mo_coeffs.shape:
                                    fix_orbital_phase(new_data.mo_coeff, old_mo_coeffs)
                        old_mo_coeffs = np.copy(new_data.mo_coeff) # for comparison of orbital phases
                        mesh1, mesh2 = mo_batch(new_data, orb_index=orb_index,spin_idx=spin, 
                                iso_level=iso_level,nx=nx,ny=ny,nz=nz,padding=padding) 
                        pl.add_mesh(mesh1, name=f"orb_pos_{i:03d}", color=color_pos)
                        pl.add_mesh(mesh2, name=f"orb_neg_{i:03d}", color=color_neg)
                    case "esp":
                        esp_mesh, active_scalar = draw_esp_molden(new_data,iso_val=iso_level,
                                        nx=nx,ny=ny,nz=nz,padding=padding)
                        mesh_args = {
                            "mesh":esp_mesh,
                            "scalars":active_scalar,
                            "cmap":cmap, 
                            "clim":[v_min, v_max], 
                            "opacity":trans,
                            "smooth_shading":True
                        }
                        pl.add_mesh(name=f"esp_{i:03d}",**mesh_args)

                    case "spin":
                        mesh1, mesh2 = draw_spin(new_data,iso_val=iso_level,nx=nx,ny=ny,nz=nz,padding=padding)
                        pl.add_mesh(mesh1, name=f"spin_pos_{i:03d}", color=color_pos)
                        pl.add_mesh(mesh2, name=f"spin_neg_{i:03d}", color=color_neg)
                    case "spin-m":
                        esp_mesh, active_scalar = draw_spin_mapped(new_data,iso_val=iso_level,nx=nx,ny=ny,
                                                        nz=nz,padding=padding)
                        mesh_args = {
                            "mesh":esp_mesh,
                            "scalars":active_scalar,
                            "cmap":cmap, 
                            "clim":[v_min, v_max], 
                            "opacity":trans,
                            "smooth_shading":True
                            }
                        pl.add_mesh(name=f"spin-m_{i:03d}", **mesh_args)  
             
                pl.export_gltf(f"{o_file}_{i:03d}.glb")
            print(f"\nDone! Multi-file GLB written to {o_file}_*.glb")
            generate_blender_script(o_file)
            #
    elif o_mode == "bld-one":
        if i_files:
            pl = pv.Plotter(off_screen=True) # one plotter for all frames
            p_bar = tqdm(i_files)
            old_mo_coeffs = None
            
            for i, name in enumerate(p_bar, start=1):
                p_bar.set_description(f"Batch-Processing {name}")
                new_data = MoleculeData.from_molden(name)
                
                # 1. MOLECULE
                comb_mesh = draw_mol(new_data.atom_points, new_data.atom_types, cpk_colors)
                # Name: mol_001, mol_002... (important for Master-Logic)
                pl.add_mesh(comb_mesh, name=f"mol_{i:03d}", label="MAT_MOLECULE", scalars="RGB", rgb=True, smooth_shading=True)

                # 2. ADDITIONAL MESHES (Orbitals, ESP, etc.)
                # each export uses the same structure: mol, orb_pos, orb_neg, surface
                # to allow for MASTER_Object assignment in Blender
                # e.g. for "mo" surface is empty
                match type:
                    case "mo":
                        # Phase-Correction for Orbitals
                        if old_mo_coeffs is not None:
                            if isinstance(new_data.mo_coeff, (list, tuple)):
                                fix_orbital_phase(new_data.mo_coeff[spin], old_mo_coeffs[spin])
                            else:
                                if new_data.mo_coeff.shape == old_mo_coeffs.shape:
                                    fix_orbital_phase(new_data.mo_coeff, old_mo_coeffs)
                        old_mo_coeffs = np.copy(new_data.mo_coeff)
                        m1, m2 = mo_batch(new_data, orb_index=orb_index, spin_idx=spin, 
                                          iso_level=iso_level, nx=nx, ny=ny, nz=nz, padding=padding)
                        pl.add_mesh(m1, name=f"orb_pos_{i:03d}",color=color_pos)
                        pl.add_mesh(m2, name=f"orb_neg_{i:03d}",color=color_neg)
                        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]), name=f"esp{i:03d}")
                    case "esp":
                        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]), name=f"orb_pos_{i:03d}",color=color_pos)
                        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]), name=f"orb_neg_{i:03d}",color=color_neg)
                        esp_mesh, active_scalar = draw_esp_molden(new_data, iso_val=iso_level,
                                                                 nx=nx, ny=ny, nz=nz, padding=padding)
                        esp_mesh.set_active_scalars(active_scalar)
                        pl.add_mesh(esp_mesh, name=f"esp_{i:03d}", scalars=active_scalar, 
                                    cmap=cmap, clim=[v_min, v_max], opacity=trans)
                    case "spin":
                        m1, m2 = draw_spin(new_data, iso_val=iso_level, nx=nx, ny=ny, nz=nz, padding=padding)
                        pl.add_mesh(m1, name=f"spin_pos_{i:03d}", color=color_pos)
                        pl.add_mesh(m2, name=f"spin_neg_{i:03d}", color=color_neg)
                        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]), name=f"esp{i:03d}")
                    case "spin-m":
                        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]), name=f"orb_pos_{i:03d}",color=color_pos)
                        pl.add_mesh(pv.PolyData([0.0,0.0,0.0]), name=f"orb_neg_{i:03d}",color=color_neg)
                        sm_mesh, active_scalar = draw_spin_mapped(new_data, iso_val=iso_level, 
                                                                 nx=nx, ny=ny, nz=nz, padding=padding)
                        sm_mesh.set_active_scalars(active_scalar)
                        pl.add_mesh(sm_mesh, name=f"spin-m_{i:03d}", scalars=active_scalar, 
                                    cmap=cmap, clim=[v_min, v_max], opacity=trans)

            # 3. SCALEBAR 
            if type in ["esp", "spin-m"]:
                cb_group = create_3d_colorbar_group(v_min, v_max, type, cmap)
                for mesh, base_name, kwargs in cb_group:
                    # Name ohne Index-Suffix, damit es statisch bleibt
                    pl.add_mesh(mesh, name=f"{base_name}_static", **kwargs)

            # 4. EXPORT
            output_path = f"{o_file}.glb"
            pl.export_gltf(output_path)
            pl.close()
            print(f"\nDone! One-file GLB written to {output_path}")

            # 5. SCRIPT GENERATION
            generate_blender_script_one(output_path) 
            
    else:
        print(f"invalid output mode: '{o_mode}', enter 'pov', 'bld' or 'bld-one' ")

if __name__ == '__main__':
    # filter for -B Flag directly from sys.argv, 
    # before it reaches the Argument-Parser 
    if "-B" in sys.argv:
        sys.argv.remove("-B")

    try:
        n_threads = get_optimal_cores()
        lib.num_threads(n_threads)
        print(f"Auto-Config: PySCF uses {n_threads} Threads.")
    except Exception as e:
        print(f"Could not set Threads automatically: {e}")

    main()