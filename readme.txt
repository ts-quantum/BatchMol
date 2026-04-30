# BatchMol

**BatchMol** is the automated high-performance engine of the suite, designed for the volumetric analysis and 
visualization of reaction trajectories. 
It processes sequences of `molden` or ´fchk´ files to track electronic properties as a reaction unfolds.

## Features
- **Multiple Analysis Modes**: Visualize Molecular Orbitals (MO), Electrostatic Potential (ESP), 
Spin Density, and Spin Mapping.
- **High-End Rendering: 
  Direct export to POV-Ray (.inc) and Blender (.glb).
* Blender Export:
  There are two options for Blender export:
  - Single File: All mesh objects in one .glb file.
  - Series: A series of .glb files (ideal for each time step 
    of an IRC trajectory).
  Additional outputs for Blender:
  - Setup Script: A '{o-file}_setup.py' is created to automate 
    processing within Blender. (see User Guide in the script)
  - Scalebar: A separate '{o-file}-scalebar.glb' is generated 
    for ESP and Spin Mapping.
* POV-Ray Integration:
  - The .inc file can be imported into POV-Ray using 'obj_name' 
    as the basename for the array of molecule objects.
  - Easy Customization: The include file contains a dedicated 
    configuration section to adapt transparency, colors, and 
    atom/bond radii.
- **Dynamic Scaling**: Automated or manual limit handling for scalebars (`v_max`).
- **Flexible Grid Control**: Full control over grid density (`n_pts`) and spatial padding.

## Command Line Interface

Usage: `python batchmol.py [FILES] [OPTIONS]`

### Arguments
- `FILES`: Input files (e.g., `scan.*.molden`). Supports wildcards.

### Options

| Option | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `--o_file` | `-o` | `output.inc/.glb` | Output file name |
| `--obj_name` | `-n` | `molecule` | Object name in output.inc |
| `--iso_level`| `-l` | `0.02` / `0.002` | Iso level (orbitals / surfaces) |
| `--n_pts` | | `50 50 50` | Number of grid points (X Y Z) |
| `--padding` | | `4` | Padding for grid calculation |
| `--type` | `-p` | `esp` | Execution mode: mo, esp, spin, spin-m |
| `--cmap` | `-c` | `rainbow` | Colormap for ESP plot |
| `--orb_index`| `-i` | `0` | Orbital index for MO plot |
| `--spin` | `-s` | `0` | Alpha (0) or beta (1) spin |
| `--v_max` | `-v` | `0.15` / `1` | Scalebar limits (v_min = -v_max) |
| `--o_mode` | `-M` | `pov` | Mode: pov, bld (Multi), bld-one |

## Project Structure

```text
.
├── batch.py               #  Main application entry point
├── requirements.txt        # Project dependencies
├── README.md               # Documentation
├── /examples               # examples with Input/Output files
│   ├── ex1                 # 1,5 H-Shift 
│   ├── ex2                 # Bromination of Propene

## Examples
### ⚠️ Configuration Note
The provided scripts (`*_split.py`, `run_*_batch.sh`, and the Blender import scripts)
contain hardcoded paths and command names. Please adjust these to match your local environment:
[_split.py and run_*_batch.sh have been created with MolAlign]

- **Quantum Chemistry**: Ensure the commands for `orca` and the conversion tool `orca_2mkl` 
(to generate `.molden` files) are correctly defined in the shell and split scripts.
- **Blender**: Update the internal paths within the Python import scripts if your project 
structure differs from the example setup.
- **File Permissions**: On Linux/macOS, remember to make the batch script executable: 
`chmod +x run_h-shift_batch.sh`.

Example 1: 1,5-H-Shift (Full Workflow Integration) - HF/6-31G [1]

  This example demonstrates the complete pipeline from a raw IRC trajectory to 
  high-end ESP and MO animations.

  1. Data Preparation (MolAlign)
  The initial IRC trajectory ts.irc_IRC_Full_trj.xyz (from an ORCA job) is processed 
  to generate a consistent trajectory and a batch-input script:
  python3 "MolAlign.py" ts.irc_IRC_Full_trj.xyz --xyz --log -f h-shift --split orca
  This generates h-shift.xyz and the split-script h-shift_split.py.

  2. Quantum Chemical Processing
  Running h-shift_split.py creates individual ORCA input files and a shell script 
  run_h-shift_batch.sh. Executing this batch job produces the .molden files required 
  for the next step.
  All intermediate files and ORCA inputs are provided in the ex1/orca directory.

  3. Volumetric Analysis & Export (BatchMol)
  BatchMol is then used to generate four different visualization outputs 
  (POV-Ray and Blender):
  Molecular Orbitals (HOMO):
  # POV-Ray include file 
  python3 batch.py ./orca/*.molden -p mo -i 18 -M pov -n homo
  (output, rendering input and final video  in ./ex1/pov-homo)
  # Blender GLB files
  python3 batch.py ./orca/*.molden -p mo -i 18 -M bld-one
  (output in ./ex1/bld-homo, final video ./ex1/bld-homo.mp4)

  4. Electrostatic Potential (ESP) Mapping:
  # POV-Ray (using HSV colormap and fixed scaling)
  python3 batch.py ./orca/*.molden -p esp -v 0.036 -M pov -c hsv -o ESP -n esp
  (output, rendering input and final video  in ./ex1/pov-esp)
  # Blender (using HSV colormap and fixed scaling)
  python3 batch.py ./orca/*.molden -p esp -v 0.036 -M bld-one -c hsv -o ESP
  (output in ./ex1/bld-esp, final video ./ex1/bld-esp.mp4)

  5. Final Rendering
  # Blender [3]: Open 'movie_template.blend' and import the respective .glb file. 
  Run the *_setup.py script to initialize the scene. Use the Trajectory_Control 
  object for positioning and scaling, and fine-tune the appearance via the 
  Dummy* objects (note: surface properties must be defined in the script before execution). 
  Finally, import the 'ESP_scalebar.glb' separately and position it in the scene 
  before rendering the video.
  # POV-Ray [4]: Use the generated .inc files with the provided video.pov and video.ini templates 
  to render the final ray-traced frames and combine them into an MP4 video (e.g., via FFmpeg [5])

Example 2: Radical Bromination of Propene (First Step) - B3LYP/def2-SVP [1]

  This example focuses on the addition of a bromine radical to propene, 
  highlighting the shift of spin density during the C-Br bond formation.

  1. Trajectory Preparation & Reversing
  Since the initial IRC might be oriented in the reverse reaction direction, we use MolAlign to 
  synchronize and flip the trajectory:

  # -r 0: Reverse the first (and only) trajectory segment
  python3 "MolAlign.py" ts1.irc_IRC_Full_trj.xyz --xyz --log -f radical --split orca -r 0

  2. Batch Input Generation (Important Step)
  ⚠️ Configuration Required: Before running the split-script, open radical_split.py and adjust 
  the following parameters:
  DFT Method: Set the appropriate functional (e.g., B3LYP def2-SVP).
  Multiplicity: Ensure it is set correctly for the radical species (e.g., 2 for doublet).
  Then, generate and execute the batch:
  python3 radical_split.py
  chmod +x run_radical_batch.sh
  ./run_radical_batch.sh  # Generates the required .molden files

  3. Volumetric Spin Mapping (BatchMol)
  This example uses the turbo colormap to visualize the spin density evolution with a fixed scale to 
  ensure a smooth animation:
  # POV-Ray include file
  python3 batch.py ./orca/*.molden -o SPIN -n spin -p spin-m -v 0.5 -M pov -c turbo
  (output, rendering input and final video  in ./ex2/pov-spin)
  # Blender GLB files (multi file)
  python3 batch.py ./orca/*.molden -o SPIN -p spin-m -v 0.5 -M bld -c turbo
  (output in ./ex1/bld-spin, final video ./ex2/bld-spin.mp4)

  4. Final Rendering
  # POV-Ray [4]:
  Use the provided video.pov and video.ini templates to render the POV-Ray frames and combine 
  them into an MP4 video (e.g., via FFmpeg [5]).
  # Blender [3]:
  Import all SPIN_*.glb files into movie_template.blend. Adjust surface settings (e.g., Alpha or Emission) 
  within SPIN_setup.py before running the script. If necessary, delete any redundant 'Renderer Nodes' 
  in the collection. Use the Trajectory_Control object to adjust the position and the Dummy_mol object 
  for molecule properties. Finally, import and position the SPIN_scalebar.glb.

[1] F. Neese, "Software update: the ORCA program system — Version 6.0", Wiley Interdiscip. Rev.: Comput. Mol. Sci., 
15, e70019 (2025). doi: 10.1002/wcms.70019.

Example 3 1-5-H-Shift (Splitting for BatchMol) B3LYP/cc-pVDZ [2]
  Pathway Description:
  This case study explores the [1,5]-sigmatropic hydrogen migration between the terminal methyl group 
  and the carbonyl oxygen of but-2-en-1-one. This classic rearrangement serves as a perfect model to 
  study the continuous transformation of electronic structures during a chemical reaction.
  Computational Workflow:
  The Intrinsic Reaction Coordinate (IRC) trajectory was high-level mapped using PSI4. To bridge the 
  gap between raw quantum chemical data and high-end visualization, MolAlign was utilized to 
  generate a dedicated post-processing script: irc_split.py.
  Key Features & Visualization:
  Automated Splitting: The script extracts every point along the IRC path, generating individual .molden 
  and .fchk files.
  BatchMol Integration: These files are ready for seamless processing with BatchMol, allowing the 
  automated generation of input data for both Blender and POV-Ray.
  Electronic Dynamics: The provided examples demonstrate how to visualize the dynamic evolution 
  of molecular properties, such as the Highest Occupied Molecular Orbital (HOMO) and the Electrostatic 
  Potential (ESP), as the reaction progresses.
  Raytracing Excellence: 
  # POV-Ray [4]:
  Pre-configured POV-Ray *.inc files for HOMO and ESP transitions are included, 
  showcasing the transformation from a reactant to a product in publication-quality renderings.
  Use the provided video.pov and video.ini templates to render the POV-Ray frames and combine 
  them into an MP4 video (e.g., via FFmpeg [5]).
  # Blender [3]:
  Open 'movie_template.blend' and import the respective .glb file. 
  Run the *_setup.py script to initialize the scene. Use the Trajectory_Control 
  object for positioning and scaling, and fine-tune the appearance via the 
  Dummy* objects (note: surface properties must be defined in the script before execution). 
  Finally, import the 'h-shift_scalebar.glb' separately and position it in the scene 
  before rendering the video.

[2] D. G. A. Smith, L. A. Burns, et al., "Psi4 1.4: Open-source software for high-throughput quantum chemistry", 
J. Chem. Phys., 152, 184108 (2020). doi: 10.1063/5.0006002.

# Rendering Software:
[3] Blender Foundation (2026). Blender (Version 5.1): Cycles Rendering Engine [Computer software]. Retrieved from blender.org
[4] POV-Ray Team (2013). Persistence of Vision Raytracer (Version 3.7) [Computer software]. GNU Affero General Public License. Retrieved from povray.org
[5] Tomar, S. (2006). Converting video formats with FFmpeg. Linux Journal, 2006(146), 10.

## Installation

1. Clone the repository
    git clone https://github.com
    cd MolVista

2. Install dependencies (using a VENV is recommended)
    pip install -r requirements.txt

    Requirements
        Python 3.x
        PyVista
        PyScf
        ....

Note: Precompiled executables for macOS, Linux arm64, and Linux x64 are available under 'releases'.

## Usage

    python3 batch.py

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.