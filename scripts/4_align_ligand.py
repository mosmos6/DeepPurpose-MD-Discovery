# scripts/4_align_ligand.py
# Author: Iori Mochizuki + patch
# Updated: 2025-10-30
# Description: Align many ligands to their docked coordinates; write per‑ligand PDB+SDF.
# Key: uses ligands.index for names; SDF is written with RDKit from the original MOL2 topology.

from pathlib import Path
import re
import numpy as np
from rdkit import Chem

from utils import extract_coordinates, kabsch, fix_pdb_element_column

INDEX_FILE = Path("ligands.index")   # written by Step 1

def _read_names_from_index() -> list[str]:
    if INDEX_FILE.exists():
        names = [ln.strip() for ln in INDEX_FILE.read_text(encoding="utf-8").splitlines()
                 if ln.strip() and not ln.strip().startswith("#")]
        # keep order, dedupe
        seen, out = set(), []
        for n in names:
            if n not in seen:
                out.append(n); seen.add(n)
        return out
    return []

def _fallback_discover_names() -> list[str]:
    # derive names from output_*.pdb (Step 3 products), e.g., output_GMP.pdb → "GMP"
    names = []
    for p in sorted(Path(".").glob("output_*.pdb")):
        m = re.match(r"output_(.+)\.pdb$", p.name)
        if m:
            names.append(m.group(1))
    # dedupe keep-order
    seen, out = set(), []
    for n in names:
        if n not in seen:
            out.append(n); seen.add(n)
    # last fallback: legacy single 'ligand'
    if not out and (Path("output.pdb").exists() or Path("ligand.pdb").exists()):
        out = ["ligand"]
    return out

def _discover_jobs():
    """
    Return list of jobs:
      {"name": N, "orig_pdb": N.pdb, "stripped_pdbqt": N.pdbqt, "docked_pdb": output_N.pdb, "template_mol2": N.mol2}
    with robust fallbacks (legacy single-ligand).
    """
    names = _read_names_from_index()
    if not names:
        names = _fallback_discover_names()
    jobs = []

    single_mode = (len(names) == 1 and names[0] == "ligand")

    for nm in names:
        orig = Path(f"{nm}.pdb")
        pq  = Path(f"{nm}.pdbqt")
        mol2 = Path(f"{nm}.mol2")
        dock = Path(f"output_{nm}.pdb")
        if not dock.exists():
            # allow "best_docked_<name>.pdb" from Step 3
            bd = Path(f"best_docked_{nm}.pdb")
            dock = bd if bd.exists() else dock
        # legacy single-ligand fallbacks
        if not orig.exists() and single_mode and Path("ligand.pdb").exists():
            orig = Path("ligand.pdb")
        if not pq.exists() and single_mode and Path("ligand.pdbqt").exists():
            pq = Path("ligand.pdbqt")
        if not mol2.exists() and single_mode and Path("ligand.mol2").exists():
            mol2 = Path("ligand.mol2")
        if not dock.exists() and single_mode and Path("output.pdb").exists():
            dock = Path("output.pdb")

        missing = [k for k, p in {"orig_pdb":orig, "stripped_pdbqt":pq, "docked_pdb":dock, "template_mol2":mol2}.items() if not p.exists()]
        if missing:
            print(f"⚠️  Skipping {nm}: missing files {missing}")
            continue

        jobs.append({"name": nm, "orig_pdb": orig, "stripped_pdbqt": pq, "docked_pdb": dock, "template_mol2": mol2})

    if not jobs:
        raise SystemExit("No ligand triplets found for alignment. Check ligands.index and Step 3 outputs (output_<name>.pdb).")
    return jobs

def _write_sdf_from_mol2_with_coords(tag: str, mol2_path: Path, coords_xyz: np.ndarray, out_sdf: Path):
    """Load RDKit Mol from MOL2 (preserves bond orders & charges), set coords (Å), write SDF."""
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    if mol is None:
        raise RuntimeError(f"[{tag}] RDKit failed to read MOL2: {mol2_path}")
    n = mol.GetNumAtoms()
    if n != coords_xyz.shape[0]:
        raise RuntimeError(f"[{tag}] Atom count mismatch: MOL2 has {n}, aligned coords have {coords_xyz.shape[0]}")
    if mol.GetNumConformers() == 0:
        conf = Chem.Conformer(n)
        mol.AddConformer(conf, assignId=True)
    conf = mol.GetConformer(0)
    for i in range(n):
        x, y, z = float(coords_xyz[i,0]), float(coords_xyz[i,1]), float(coords_xyz[i,2])
        conf.SetAtomPosition(i, (x, y, z))
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol); w.close()

def _align_one(job):
    nm = job["name"]
    docked_coords, _   = extract_coordinates(str(job["docked_pdb"]))
    full_coords, full_atoms = extract_coordinates(str(job["orig_pdb"]))
    stripped_coords, _ = extract_coordinates(str(job["stripped_pdbqt"]))

    print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms for {nm}")

    # stabilize pairing: sort by radial distance to centroid
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned = (full_coords - c_stripped) @ R.T + c_docked  # Å

    aligned_pdb = Path(f"aligned_{nm}_fixed.pdb")
    fixed_pdb   = Path(f"fixed_{nm}.pdb")
    out_sdf     = Path(f"{nm}.sdf")

    # PDB for QC
    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")

    fix_pdb_element_column(str(aligned_pdb), str(fixed_pdb))
    _write_sdf_from_mol2_with_coords(nm, job["template_mol2"], aligned, out_sdf)
    print(f"✅ {nm}: aligned → {aligned_pdb} | fixed → {fixed_pdb} | SDF → {out_sdf}")

def main():
    jobs = _discover_jobs()
    print(f"\n[4] Aligning {len(jobs)} ligand(s): {', '.join(j['name'] for j in jobs)}\n")
    for j in jobs:
        _align_one(j)

    # index of SDFs for Step 5 (and write a friendly summary)
    names = [j["name"] for j in jobs]
    Path("ligands_sdf.index").write_text("\n".join(names) + "\n", encoding="utf-8")
    print("\nSummary:")
    print("  • ligands_sdf.index  (one name per line, SDF = <name>.sdf)")
    for n in names:
        print(f"    - {n}.sdf")

if __name__ == "__main__":
    main()
