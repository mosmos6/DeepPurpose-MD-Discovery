### 3_docking_vina.py
# Author: Iori Mochizuki
# Updated: 2025-09-17
# Description: Run docking using AutoDock Vina and extract best pose

"""
Step 3: Docking with AutoDock Vina
- Optionally accept a manual box center via --center_x/--center_y/--center_z
- Otherwise:
  * Calculate centroid (box center) from receptor structure
  * Use RNA residue-range centroid if requested
  * Or use HETATM centroid (if present and --use-residue-centroid)
  * Or fallback to all-ATOM centroid
- Run AutoDock Vina and extract best pose (MODEL 1)
"""

import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--use-residue-centroid", action="store_true",
                    help="Use residue-based centroid from HETATM entries if available.")
parser.add_argument("--rna", action="store_true",
                    help="Enable RNA mode (calculate centroid on nucleotides)")
parser.add_argument("--res-range", type=str, default=None,
                    help="Residue range (e.g. 13-41) to calculate centroid (for RNA)")
# NEW: manual center override
parser.add_argument("--center_x", type=float, default=None, help="Manual docking box center X (Å)")
parser.add_argument("--center_y", type=float, default=None, help="Manual docking box center Y (Å)")
parser.add_argument("--center_z", type=float, default=None, help="Manual docking box center Z (Å)")

args = parser.parse_args()

receptor_file = "receptor_for_centroid.pdb"

def parse_res_range(res_range_str):
    if res_range_str:
        for dash in ['–', '—', '−', '‒', '―']:
            res_range_str = res_range_str.replace(dash, '-')
        start, end = [int(x) for x in res_range_str.split('-')]
        return start, end
    return None, None

# --------------------------------------------
# Priority 1: manual override if all provided
# --------------------------------------------
manual_center = (
    args.center_x is not None and
    args.center_y is not None and
    args.center_z is not None
)

if manual_center:
    x_center = round(float(args.center_x), 3)
    y_center = round(float(args.center_y), 3)
    z_center = round(float(args.center_z), 3)
    print(f"📍 Docking box center (manual override): x={x_center}, y={y_center}, z={z_center}")

# --------------------------------------------
# Priority 2: RNA residue-range centroid
# --------------------------------------------
elif args.rna and args.res_range:
    print(f"🔎 RNA MODE: Calculating centroid for residues {args.res_range} ...")
    start, end = parse_res_range(args.res_range)
    x_sum = y_sum = z_sum = 0.0
    atom_count = 0
    with open(receptor_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                resid = int(line[22:26].strip())
                if start <= resid <= end:
                    x_sum += float(line[30:38])
                    y_sum += float(line[38:46])
                    z_sum += float(line[46:54])
                    atom_count += 1
    if atom_count == 0:
        raise ValueError("No atoms found in the specified residue range!")
    x_center = round(x_sum / atom_count, 3)
    y_center = round(y_sum / atom_count, 3)
    z_center = round(z_sum / atom_count, 3)
    print(f"📍 RNA box center (residues {start}-{end}): x={x_center}, y={y_center}, z={z_center}")

# --------------------------------------------
# Priority 3: HETATM centroid (residue-based)
# --------------------------------------------
elif args.use_residue_centroid:
    print("🔎 Calculating centroid from HETATM residues...")
    residue_names = set()
    with open(receptor_file, 'r') as f:
        for line in f:
            if line.startswith("HETATM"):
                residue_names.add(line[17:20].strip())
    target_residues = list(residue_names)
    print(f"Detected target residues: {target_residues}")

    x_sum = y_sum = z_sum = 0.0
    atom_count = 0
    with open(receptor_file, 'r') as f:
        for line in f:
            if line.startswith("HETATM") and any(res in line for res in target_residues):
                x_sum += float(line[30:38])
                y_sum += float(line[38:46])
                z_sum += float(line[46:54])
                atom_count += 1

    if atom_count == 0:
        print("⚠️ No HETATM residues found — falling back to ATOM centroid.")
        args.use_residue_centroid = False
    else:
        x_center = round(x_sum / atom_count, 3)
        y_center = round(y_sum / atom_count, 3)
        z_center = round(z_sum / atom_count, 3)
        print(f"📍 Docking box center (residue centroid): x={x_center}, y={y_center}, z={z_center}")

# --------------------------------------------
# Priority 4: fallback ATOM centroid
# --------------------------------------------
if not manual_center and not (args.rna and args.res_range) and not args.use_residue_centroid:
    x_sum = y_sum = z_sum = 0.0
    atom_count = 0
    with open(receptor_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                x_sum += float(line[30:38])
                y_sum += float(line[38:46])
                z_sum += float(line[46:54])
                atom_count += 1
    x_center = round(x_sum / atom_count, 3)
    y_center = round(y_sum / atom_count, 3)
    z_center = round(z_sum / atom_count, 3)
    print(f"📍 Docking box center (all ATOMs): x={x_center}, y={y_center}, z={z_center}")

# --------------------------------------------
# Run AutoDock Vina
# --------------------------------------------
vina_command = f"""
./vina_1.2.5_linux_x86_64 \
--receptor receptor.pdbqt \
--ligand ligand.pdbqt \
--out output.pdbqt \
--center_x {x_center} \
--center_y {y_center} \
--center_z {z_center} \
--size_x 30 --size_y 30 --size_z 30 \
--seed 12345 --exhaustiveness 20
"""
print("Running AutoDock Vina...")
subprocess.run(vina_command, shell=True, check=True)

# Extract MODEL 1 as best pose
with open("output.pdbqt", "r") as infile, open("best_docked_ligand.pdb", "w") as outfile:
    inside_model = False
    for line in infile:
        if line.startswith("MODEL 2"):
            break
        if line.startswith("MODEL 1"):
            inside_model = True
        if inside_model and not line.startswith("MODEL"):
            outfile.write(line)

print("✅ Best docking pose extracted → best_docked_ligand.pdb")

# Clean docking-specific columns and write final PDB
with open("best_docked_ligand.pdb", "r") as infile, open("output.pdb", "w") as outfile:
    for line in infile:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            outfile.write(line)

print("✅ Cleaned ligand pose written to output.pdb")
