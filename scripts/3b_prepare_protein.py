# 3b_prepare_protein.py
# Author: Iori Mochizuki
# Updated: 2025-07-25
# Description: Clean and standardize receptor (protein or RNA) using PDBFixer. Remove terminal phosphates for RNA.

import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller, ForceField
import os

parser = argparse.ArgumentParser(description="Clean/standardize protein or RNA receptor for MD.")
parser.add_argument("--rna", action="store_true", help="Enable RNA-aware mode (SimRNA input, special patching)")
parser.add_argument("--input-pdb", type=str, default="receptor_clean.pdb", help="Input PDB (default: receptor_clean.pdb)")
parser.add_argument("--output-pdb", type=str, default="receptor_cleaned.pdb", help="Output PDB (default: receptor_cleaned.pdb)")
args = parser.parse_args()

def get_residue_range(pdb_path):
    first_res = None
    last_res = None
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_id = int(line[22:26].strip())
                if first_res is None:
                    first_res = res_id
                last_res = res_id  # keeps updating
    return f"{first_res:4d}", f"{last_res:4d}"  # Properly padded

if args.rna:
    print("🧬 [RNA MODE] Running PDBFixer gently (clear missing residues)...")
    fixer = PDBFixer(filename=args.input_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # <-- Key: do NOT add residues!
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    # Remove heterogens except water
    fixer.removeHeterogens(keepWater=True)
    # Write temp
    with open("receptor_fixed.pdb", "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print("✅ PDBFixer complete: receptor_fixed.pdb")

    # Patch 5' terminal phosphates
    print("🩺 Removing 5′ terminal phosphates from first residue...")
    input_pdb = "receptor_fixed.pdb"
    output_pdb = "receptor_patched.pdb"
    first_res_id, last_res_id = get_residue_range(input_pdb)
    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                resnum = line[22:26]
                atom_name = line[12:16].strip()
                # Remove 5′ terminal phosphates from first residue
                if resnum == first_res_id and atom_name in ["P", "OP1", "OP2"]:
                    continue
            fout.write(line)
    print(f"✅ Terminal patching complete → {output_pdb}")

    # (Optional: hydrogenate with RNA forcefield to fix atom names, ensure all H's present)
    print("💧 Adding hydrogens with RNA forcefield...")
    patched = PDBFile(output_pdb)
    modeller = Modeller(patched.topology, patched.positions)
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
    modeller.addHydrogens(forcefield)
    with open(args.output_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"🎉 Final RNA receptor ready: {args.output_pdb}")

else:
    # --- Protein: Standard PDBFixer workflow ---
    fixer = PDBFixer(filename=args.input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    fixer.removeHeterogens(keepWater=True)
    with open(args.output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"✅ Protein receptor ready: {args.output_pdb}")
