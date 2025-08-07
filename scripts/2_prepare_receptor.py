# 2_prepare_receptor.py
# Author: Iori Mochizuki
# Updated: 2025-07-25
# Description: Prepare protein (by PDB ID) or RNA (SimRNA PDB) for docking and MD

import os, sys, subprocess, argparse

parser = argparse.ArgumentParser(description="Prepare receptor (protein or RNA) for docking and MD.")
parser.add_argument("--pdb-id", type=str, default=None, help="PDB ID to download from RCSB (protein mode)")
parser.add_argument("--rna", action="store_true", help="Enable RNA mode: use uploaded PDB file (SimRNA output)")
parser.add_argument("--input-pdb", type=str, default="receptor.pdb", help="Input RNA PDB file if --rna (default: receptor.pdb)")
parser.add_argument("--skip-fix", action="store_true", help="Skip PDBFixer step (protein mode only)")
parser.add_argument("--strict-protein", action="store_true", help="Strip non-protein residues using MDTraj (protein mode only)")
args = parser.parse_args()

receptor_clean = "receptor_clean.pdb"
receptor_for_centroid = "receptor_for_centroid.pdb"

if args.rna:
    # ---- RNA MODE ----
    print("🧬 [RNA MODE] Using uploaded PDB:", args.input_pdb)
    # Copy the input RNA PDB to standard names
    subprocess.run(["cp", args.input_pdb, receptor_clean])
    subprocess.run(["cp", args.input_pdb, receptor_for_centroid])
    receptor_for_meeko = receptor_clean
else:
    # ---- PROTEIN MODE ----
    if not args.pdb_id:
        print("❌ Error: --pdb-id must be specified unless --rna is set.")
        sys.exit(1)
    receptor_raw = "receptor.pdb"
    receptor_fixed = "receptor_fixed.pdb"
    print(f"📦 Downloading {args.pdb_id} from RCSB...")
    subprocess.run(["wget", f"https://files.rcsb.org/download/{args.pdb_id}.pdb", "-O", receptor_raw])

    # Remove water
    print("🧹 Removing water molecules...")
    subprocess.run(["obabel", receptor_raw, "-O", receptor_clean, "--delete", "HOH"])
    subprocess.run(["obabel", receptor_raw, "-O", receptor_for_centroid, "--delete", "HOH"])

    # Optional MDTraj cleaning
    if args.strict_protein:
        print("🔬 Removing non-protein residues using MDTraj...")
        import mdtraj as md
        traj = md.load_pdb(receptor_clean)
        filtered = traj.atom_slice(traj.top.select("protein"))
        filtered.save_pdb(receptor_clean)
        print("✅ Saved as receptor_clean.pdb")

    # Optional PDBFixer
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

# --- Meeko: Always run for both protein and RNA ---
print("🛠️ Running Meeko on:", receptor_for_meeko)
meeko_command = (
    "python3 /usr/local/envs/deeppurpose-md-env/lib/python3.11/site-packages/meeko/cli/mk_prepare_receptor.py "
    f"-i {receptor_for_meeko} -o receptor -p -j -v "
    "--box_size 20 20 20 --box_center 0 0 0 --allow_bad_res"
)
os.system(meeko_command)
print("✅ Receptor PDBQT prepared.")
