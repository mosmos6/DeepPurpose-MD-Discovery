"""
Script Name: 1_prepare_ligand.py
Description: Generates 3D ligand structure and formats it for docking.
Author: Iori Mochizuki
Date: 2025-05-13
"""

import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

# --- SETTINGS ---
smiles = sys.argv[1]  
ligand_name = "ligand"

# --- 1. Generate RDKit molecule from SMILES ---
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(mol)

# --- 2. Write PDB file ---
pdb_path = f"{ligand_name}.pdb"
with open(pdb_path, 'w') as f:
    f.write(Chem.MolToPDBBlock(mol))
print(f"✅ PDB file written: {pdb_path}")

# --- 3. Convert to mol2 using Open Babel ---
mol2_path = f"{ligand_name}.mol2"
subprocess.run(["obabel", pdb_path, "-O", mol2_path])
print(f"✅ mol2 file created: {mol2_path}")

# --- 4. Convert to PDBQT for Vina docking ---
pdbqt_path = f"{ligand_name}.pdbqt"
subprocess.run(["obabel", mol2_path, "-O", pdbqt_path])
print(f"✅ PDBQT file created: {pdbqt_path}")
