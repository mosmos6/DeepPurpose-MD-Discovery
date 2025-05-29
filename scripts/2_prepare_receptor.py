# 2_prepare_receptor.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Fetch and prepare receptor PDB file for docking and MD

"""
Step 2: Receptor Preparation
- Download PDB structure from RCSB
- Remove water molecules
- Optionally skip PDBFixer for Meeko-sensitive structures
- Output ready-to-use receptor files for docking
"""


import os, sys, subprocess, argparse

parser = argparse.ArgumentParser()
parser.add_argument("pdb_id", type=str, help="PDB ID to download")
parser.add_argument("--skip-fix", action="store_true", help="Skip PDBFixer step")
parser.add_argument("--strict-protein", action="store_true", help="Strip non-protein residues using MDTraj")
args = parser.parse_args()

PDB_ID = args.pdb_id
receptor_raw = "receptor.pdb"
receptor_clean = "receptor_clean.pdb"
receptor_fixed = "receptor_fixed.pdb"
receptor_for_centroid = "receptor_for_centroid.pdb"

# --- Download ---
print(f"📦 Downloading {PDB_ID}...")
subprocess.run(["wget", f"https://files.rcsb.org/download/{PDB_ID}.pdb", "-O", receptor_raw])

# --- Remove water ---
print("🧹 Removing water molecules...")
subprocess.run(["obabel", receptor_raw, "-O", receptor_clean, "--delete", "HOH"])
subprocess.run(["obabel", receptor_raw, "-O", receptor_for_centroid, "--delete", "HOH"])

# --- Optional MDTraj cleaning ---
if args.strict_protein:
    print("🔬 Removing non-protein residues using MDTraj...")
    import mdtraj as md
    traj = md.load_pdb(receptor_clean)
    filtered = traj.atom_slice(traj.top.select("protein"))
    filtered.save_pdb("receptor_clean.pdb")  # overwrite
    print("✅ Saved as receptor_clean.pdb")

# --- PDBFixer ---
if args.skip_fix:
    print("⚠️ Skipping PDBFixer cleanup (as requested via --skip-fix).")
    receptor_for_meeko = receptor_clean
else:
    print("🧬 Running PDBFixer...")
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=receptor_clean)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    with open(receptor_fixed, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"✅ PDBFixer output → {receptor_fixed}")
    receptor_for_meeko = receptor_fixed

# --- Meeko ---
print("🛠️ Running Meeko on:", receptor_for_meeko)
meeko_command = (
    f"python3 /usr/local/envs/deeppurpose-md-env/lib/python3.11/site-packages/meeko/cli/mk_prepare_receptor.py "
    f"-i {receptor_for_meeko} -o receptor -p -j -v "
    "--box_size 20 20 20 --box_center 0 0 0 --allow_bad_res"
)
os.system(meeko_command)
print("✅ Receptor PDBQT prepared.")
