"""
Script Name: 1_prepare_ligand.py
Description: Generate 3D ligands from one or many SMILES and format for docking.
Author: Iori Mochizuki + patch
Updated: 2025-10-28

Accepted inputs (all via sys.argv):
  1) Single SMILES string
     python scripts/1_prepare_ligand.py "CCO"

  2) JSON/Python dict mapping names->smiles
     python scripts/1_prepare_ligand.py '{"aspirin":"CC(=O)OC1=CC=CC=C1C(=O)O","caffeine":"Cn1cnc2n(C)c(=O)n(C)c(=O)c12"}'

  3) JSON/Python list of SMILES (auto-named ligand_1, ligand_2, ...)
     python scripts/1_prepare_ligand.py '["CCO","c1ccccc1"]'

  4) Comma-separated string (auto-named ligand_1, ligand_2, ...)
     python scripts/1_prepare_ligand.py "CCO,c1ccccc1"

  5) @file (one SMILES per line, optional 2nd column = name)
     python scripts/1_prepare_ligand.py @smiles.txt
"""

import sys, json, ast, re, subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# ---------- parsing helpers ----------
NAME_SAFE = re.compile(r"[^A-Za-z0-9_.-]+")

def _safe_name(s: str) -> str:
    s2 = NAME_SAFE.sub("_", s.strip())
    return s2 or "ligand"

def _parse_payload(tokens):
    """
    tokens: list of argv[1:] strings.
    Returns list of (name, smiles) pairs.
    """
    if not tokens:
        sys.exit("No SMILES given. See header for usage.")

    # Join tokens to let users pass dict/list split across args
    payload = " ".join(tokens).strip()

    # @file support
    if payload.startswith("@"):
        rows = []
        path = Path(payload[1:])
        for ln in path.read_text(encoding="utf-8").splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.replace("\t", " ").split()
            smi = parts[0]
            nm  = parts[1] if len(parts) > 1 else None
            rows.append((_safe_name(nm) if nm else None, smi))
        # name assignment below
        return _assign_names(rows)

    # Try JSON → literal → fallbacks
    obj = None
    try:
        obj = json.loads(payload)
    except Exception:
        try:
            obj = ast.literal_eval(payload)
        except Exception:
            obj = None

    if isinstance(obj, dict):
        rows = [(_safe_name(k), str(v)) for k, v in obj.items()]
        return rows

    if isinstance(obj, (list, tuple)):
        rows = [(None, str(s)) for s in obj]
        return _assign_names(rows)

    # Comma-separated “a,b,c”
    if "," in payload and " " not in payload:
        rows = [(None, p.strip()) for p in payload.split(",") if p.strip()]
        return _assign_names(rows)

    # Single SMILES string → legacy 'ligand'
    return [("ligand", payload)]

def _assign_names(pairs):
    """Fill None names as ligand_1, ligand_2,..."""
    out = []
    k = 1
    for nm, smi in pairs:
        if nm is None:
            nm = f"ligand_{k}"
            k += 1
        out.append((_safe_name(nm), smi))
    return out

# ---------- chemistry helpers ----------
def _mol_from_smiles_3d(smi: str, seed: int = 1337):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    # More robust embedding for tricky ligands
    status = AllChem.EmbedMolecule(m, params)
    if status != 0:
        # fallback: random seedless try
        AllChem.EmbedMolecule(m, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(m, maxIters=800)
    return m

def _write_one(name: str, mol):
    pdb = Path(f"{name}.pdb")
    pdb.write_text(Chem.MolToPDBBlock(mol), encoding="utf-8")
    print(f"✅ PDB written: {pdb}")

    # Open Babel conversions
    mol2 = Path(f"{name}.mol2")
    pdbqt = Path(f"{name}.pdbqt")
    subprocess.run(["obabel", str(pdb), "-O", str(mol2)], check=True)
    print(f"✅ MOL2 written: {mol2}")
    subprocess.run(["obabel", str(mol2), "-O", str(pdbqt)], check=True)
    print(f"✅ PDBQT written: {pdbqt}")

# ---------- main ----------
def main():
    items = _parse_payload(sys.argv[1:])  # list[(name, smiles)]
    if not items:
        sys.exit("No valid ligands after parsing.")

    # ensure unique names
    seen = {}
    unique_items = []
    for nm, smi in items:
        base = nm
        while nm in seen:
            seen[nm] += 1
            nm = f"{base}_{seen[nm]}"
        seen.setdefault(nm, 0)
        unique_items.append((nm, smi))

    print(f"\n[1] Preparing {len(unique_items)} ligand(s):")
    for i, (nm, smi) in enumerate(unique_items, start=1):
        print(f"  - ({i}) {nm}: {smi}")

    # build each
    for i, (nm, smi) in enumerate(unique_items, start=1):
        mol = _mol_from_smiles_3d(smi, seed=1337 + i)
        _write_one(nm, mol)

    # Back-compat: if exactly one and its name != 'ligand', also write legacy basename
    if len(unique_items) == 1 and unique_items[0][0] != "ligand":
        base = unique_items[0][0]
        for ext in (".pdb", ".mol2", ".pdbqt"):
            src = Path(f"{base}{ext}")
            if src.exists():
                dst = Path(f"ligand{ext}")
                dst.write_bytes(src.read_bytes())
                print(f"↪️  Also wrote legacy {dst} (copy of {src})")

    # Write index files for later steps
    Path("ligands.json").write_text(
        json.dumps({nm: smi for nm, smi in unique_items}, indent=2),
        encoding="utf-8"
    )
    Path("ligands.index").write_text(
        "\n".join(nm for nm, _ in unique_items) + "\n",
        encoding="utf-8"
    )
    print("\nSummary:")
    print("  • ligands.json   (name → SMILES)")
    print("  • ligands.index  (one name per line, in the creation order)")

if __name__ == "__main__":
    main()
