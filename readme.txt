# BatchMol

**BatchMol** is the automated high-performance engine of the suite, designed for the volumetric analysis and visualization of reaction trajectories. It processes sequences of `molden` files to track electronic properties as a reaction unfolds.

## Features
- **Multiple Analysis Modes**: Visualize Molecular Orbitals (MO), Electrostatic Potential (ESP), Spin Density, and Spin Mapping.
- **High-End Rendering**: Direct export to **POV-Ray (.inc)** and **Blender (.glb)**.
- **Dynamic Scaling**: Automated or manual limit handling for scalebars (`v_max`).
- **Flexible Grid Control**: Full control over grid density (`n_pts`) and spatial padding.

## Command Line Interface

Usage: `python batchmol.py [FILES] [OPTIONS]`

### Arguments
- `FILES`: Input files (e.g., `scan.*.molden`). Supports wildcards.

### Options

| Option | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `--type` | `-p` | `esp` | Mode: `mo`, `esp`, `spin`, `spin-m` (Spin Mapping) |
| `--o_mode` | `-m` | `pov` | Output: `pov` (POV-Ray) or `bld` (Blender) |
| `--iso_level`| `-l` | `0.02` | Iso-surface value |
| `--n_pts` | | `50 50 50` | Number of grid points (X Y Z) |
| `--padding` | | `4` | Grid calculation padding |
| `--cmap` | `-c` | `rainbow` | Colormap (e.g., `viridis`, `coolwarm`, `RdBu`) |
| `--orb_index`| `-i` | `0` | Orbital index for MO plots |
| `--v_max` | `-v` | *auto* | Scalebar limits (± value) |
| `--o_file` | `-o` | `output` | Base name for output files |

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
The provided scripts (`h-shift_split.py`, `run_h-shift_batch.sh`, and the Blender import scripts)
contain hardcoded paths and command names. Please adjust these to match your local environment:

- **Quantum Chemistry**: Ensure the commands for `orca` and the conversion tool `orca_2mkl` 
(to generate `.molden` files) are correctly defined in the shell and split scripts.
- **Blender**: Update the internal paths within the Python import scripts if your project 
structure differs from the example setup.
- **File Permissions**: On Linux/macOS, remember to make the batch script executable: 
`chmod +x run_h-shift_batch.sh`.

Example 1: 1,5-H-Shift (Full Workflow Integration)

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
  python3 batch.py ./orca/*.molden -p mo -i 18 -m pov -n homo
  # Blender GLB files
  python3 batch.py ./orca/*.molden -p mo -i 18 -m bld
  Electrostatic Potential (ESP) Mapping:

  # POV-Ray (using HSV colormap and fixed scaling)
  python3 batch.py ./orca/*.molden -p esp -v 0.036 -m pov -c hsv -o ESP -n esp
  # Blender (using HSV colormap and fixed scaling)
  python3 batch.py ./orca/*.molden -p esp -v 0.036 -m bld -c hsv -o ESP

  4. Final Rendering
  Blender: Import the generated .glb files using the provided import_and_animate.py script.
  POV-Ray: Use the generated .inc files with the provided video.pov and video.ini templates 
  to render the final ray-traced frames.

Example 2: Radical Bromination of Propene (First Step)

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
  DFT Method: Set the appropriate functional (e.g., UM06-2X).
  Multiplicity: Ensure it is set correctly for the radical species (e.g., 2 for doublet).
  Then, generate and execute the batch:
  python3 radical_split.py
  chmod +x run_radical_batch.sh
  ./run_radical_batch.sh  # Generates the required .molden files

  3. Volumetric Spin Mapping (BatchMol)
  We use the turbo colormap to visualize the spin density evolution with a fixed scale to 
  ensure a smooth animation:
  # POV-Ray include file
  python3 batch.py ./orca/*.molden -o SPIN -n spin -p spin-m -v 0.5 -m pov -c turbo

  # Blender GLB files
  python3 batch.py ./orca/*.molden -o SPIN -p spin-m -v 0.5 -m bld -c turbo

  4. Final Rendering
  Use the provided video.pov and video.ini templates to render the POV-Ray frames and combine 
  them into an MP4 video (e.g., via FFmpeg).

## Installation

1. Clone the repository
    git clone https://github.com
    cd MolVista

2. Install dependencies
    pip install -r requirements.txt

    Requirements
        Python 3.x
        PyVista
        PyScf
        ....

## Usage

    python3 batch.py

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.