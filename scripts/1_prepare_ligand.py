"""
Script Name: 1_prepare_ligand.py
Description: Generate 3D ligands from one or many SMILES and format for docking.
Author: Iori Mochizuki + patch
Updated: 2025-10-22

Usage examples
--------------
# single (legacy behavior preserved; writes ligand.{pdb,mol2,pdbqt})
python scripts/1_prepare_ligand.py "CCO"

# multiple as separate args (writes ligand_1.*, ligand_2.*, ...)
python scripts/1_prepare_ligand.py "CCO" "c1ccccc1"

# comma-separated in one arg
python scripts/1_prepare_ligand.py "CCO,c1ccccc1"

# from file (one SMILES per line, optional 2nd column = name)
python scripts/1_prepare_ligand.py @smiles.txt
"""

import argparse, sys, subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

def _parse_smiles_args(args_list):
    """Return list of (name, smiles). Accept @file, comma-separated, or many args."""
    out = []
    tokens = []
    for tok in args_list:
        tok = tok.strip()
        if not tok:
            continue
        if tok.startswith("@"):
            # read lines:  SMILES [name]
            for ln in Path(tok[1:]).read_text(encoding="utf-8").splitlines():
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.replace("\t", " ").split()
                smi = parts[0]
                nm = parts[1] if len(parts) > 1 else None
                tokens.append((nm, smi))
        else:
            # allow comma-separated blob
            for sub in tok.split(","):
                sub = sub.strip()
                if sub:
                    tokens.append((None, sub))
    # assign default names if absent
    if len(tokens) == 1:
        nm = tokens[0][0] or "ligand"
        out.append((nm, tokens[0][1]))
    else:
        for i, (nm, smi) in enumerate(tokens, start=1):
            out.append((nm or f"ligand_{i}", smi))
    return out

def _rdkit_3d_from_smiles(smi, seed=1337):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    AllChem.EmbedMolecule(m, params)
    AllChem.UFFOptimizeMolecule(m, maxIters=500)
    return m

def _write_files(name, mol):
    pdb = Path(f"{name}.pdb")
    pdb.write_text(Chem.MolToPDBBlock(mol), encoding="utf-8")
    print(f"✅ PDB written: {pdb}")

    # Open Babel conversions (mol2 → pdbqt)
    mol2 = Path(f"{name}.mol2")
    pdbqt = Path(f"{name}.pdbqt")
    subprocess.run(["obabel", str(pdb), "-O", str(mol2)], check=True)
    print(f"✅ MOL2 written: {mol2}")
    subprocess.run(["obabel", str(mol2), "-O", str(pdbqt)], check=True)
    print(f"✅ PDBQT written: {pdbqt}")

def main():
    ap = argparse.ArgumentParser(description="Prepare one or many ligands from SMILES.")
    ap.add_argument("smiles", nargs="+",
                    help="SMILES strings; comma‑separated OK. Or @file with one SMILES per line.")
    args = ap.parse_args()

    items = _parse_smiles_args(args.smiles)
    if not items:
        sys.exit("No valid SMILES given.")

    # Generate each ligand
    for idx, (name, smi) in enumerate(items, start=1):
        print(f"--- [{idx}/{len(items)}] {name}  SMILES={smi}")
        mol = _rdkit_3d_from_smiles(smi, seed=1337 + idx)
        _write_files(name, mol)

    # Backward compatibility: if exactly one ligand named not 'ligand',
    # also write the legacy filenames so Step 3 & 4 continue to work unchanged.
    if len(items) == 1 and items[0][0] != "ligand":
        base = items[0][0]
        for ext in (".pdb", ".mol2", ".pdbqt"):
            src = Path(f"{base}{ext}")
            if src.exists():
                dst = Path(f"ligand{ext}")
                dst.write_bytes(src.read_bytes())
                print(f"↪️  Also wrote legacy {dst} (copy of {src})")

if __name__ == "__main__":
    main()
