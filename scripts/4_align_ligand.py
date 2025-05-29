# 4_align_ligand.py
# Author: Iori Mochizuki
# Updated: 2025-05-13
# Description: Align ligand to docked coordinates using Kabsch algorithm and output SDF for OpenMM.

import numpy as np
import subprocess
from utils import extract_coordinates, kabsch, fix_pdb_element_column

# === File paths ===
docked_pdb = "output.pdb"                # Output from AutoDock Vina
original_pdb = "ligand.pdb"              # RDKit-style full ligand
stripped_pdbqt = "ligand.pdbqt"          # Ligand used for docking (used for atom matching)
aligned_pdb = "aligned_ligand_fixed.pdb" # Output before cleaning
final_pdb = "fixed_ligand.pdb"           # Final cleaned PDB
final_sdf = "ligand.sdf"                 # Output for OpenFF (used in MD)

# === Step 1: Load coordinates ===
docked_coords, docked_atoms = extract_coordinates(docked_pdb)
full_coords, full_atoms = extract_coordinates(original_pdb)
stripped_coords, _ = extract_coordinates(stripped_pdbqt)

print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms.")

# === Step 2: Sort atoms by radial distance to centroid (for Kabsch stability) ===
d1 = np.linalg.norm(docked_coords - np.mean(docked_coords, axis=0), axis=1)
d2 = np.linalg.norm(stripped_coords - np.mean(stripped_coords, axis=0), axis=1)
sorted_docked = docked_coords[np.argsort(d1)]
sorted_stripped = stripped_coords[np.argsort(d2)]

# === Step 3: Apply Kabsch alignment ===
R, centroid_stripped, centroid_docked = kabsch(sorted_stripped, sorted_docked)
aligned_coords = (full_coords - centroid_stripped) @ R.T + centroid_docked

# === Step 4: Write aligned structure ===
with open(aligned_pdb, "w") as f:
    for i, (_, _, _, _, line) in enumerate(full_atoms):
        x, y, z = aligned_coords[i]
        aligned_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        f.write(aligned_line)

print(f"✅ Aligned ligand written to: {aligned_pdb}")

# === Step 5: Fix element column for OpenMM ===
fix_pdb_element_column(aligned_pdb, final_pdb)

# === Step 6: Convert to SDF (for OpenFF ligand processing) ===
subprocess.run(["obabel", final_pdb, "-O", final_sdf])
print(f"✅ Ligand exported to {final_sdf} for forcefield assignment.")
