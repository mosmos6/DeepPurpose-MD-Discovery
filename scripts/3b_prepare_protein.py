# 3b_prepare_protein.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Clean and standardize receptor using PDBFixer

import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile

"""
Step 3b: Receptor Preprocessing
- Fix missing atoms and residues
- Replace nonstandard residues
- Add hydrogens
- Remove non-protein heterogens (ions, ligands)
"""

# Input and output filenames
INPUT_PDB = "receptor_clean.pdb"      # PDB after HOH removal
OUTPUT_PDB = "receptor_cleaned.pdb"   # Final cleaned PDB for docking & MD

# Initialize fixer
fixer = PDBFixer(filename=INPUT_PDB)

# Repair structure
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

# Remove heterogens except water
fixer.removeHeterogens(keepWater=True)

# Write cleaned structure
with open(OUTPUT_PDB, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print(f"✅ Receptor preprocessing complete: {OUTPUT_PDB}")
