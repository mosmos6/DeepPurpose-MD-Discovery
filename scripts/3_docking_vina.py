### 3_docking_vina.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Run docking using AutoDock Vina and extract best pose

"""
Step 3: Docking with AutoDock Vina
- Calculate centroid (box center) from receptor structure
- Run AutoDock Vina with specified box center and size
- Extract the best scoring pose (MODEL 1)
- Save to clean PDB file for MD alignment
"""

import subprocess
import argparse

# CLI flag: --use-residue-centroid
parser = argparse.ArgumentParser()
parser.add_argument("--use-residue-centroid", action="store_true", help="Use residue-based centroid from HETATM entries if available.")
args = parser.parse_args()

# File path
receptor_file = "receptor_for_centroid.pdb"

# === Optional: Residue-based centroid calculation ===
if args.use_residue_centroid:
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

# === Fallback: ATOM centroid (default) ===
if not args.use_residue_centroid:
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
    print(f"📍 Docking box center (default centroid): x={x_center}, y={y_center}, z={z_center}")


# Step 2: Run AutoDock Vina
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

# Step 3: Extract MODEL 1 as best pose
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

# Step 4: Clean up docking-specific columns and write final PDB
with open("best_docked_ligand.pdb", "r") as infile, open("output.pdb", "w") as outfile:
    for line in infile:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            outfile.write(line)

print("✅ Cleaned ligand pose written to output.pdb")
