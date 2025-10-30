# Author: Iori Mochizuki + patch
# Updated: 2025-10-30
# Description: Align many ligands to their docked coordinates; write per‑ligand PDB+SDF.
# Key change: SDF is written with Open Babel from the fixed aligned PDB (robust),
#             with an RDKit fallback (sanitize=False) only if OBabel fails.

from pathlib import Path
import re, subprocess
import numpy as np

from utils import extract_coordinates, kabsch, fix_pdb_element_column

INDEX_FILE = Path("ligands.index")  # written by Step 1

def _read_names_from_index() -> list[str]:
    if INDEX_FILE.exists():
        names = [ln.strip() for ln in INDEX_FILE.read_text(encoding="utf-8").splitlines()
                 if ln.strip() and not ln.strip().startswith("#")]
        seen, out = set(), []
        for n in names:
            if n not in seen:
                out.append(n); seen.add(n)
        return out
    return []

def _fallback_discover_names() -> list[str]:
    # derive names from Step 3 products: output_<name>.pdb → "<name>"
    names = []
    for p in sorted(Path(".").glob("output_*.pdb")):
        m = re.match(r"output_(.+)\.pdb$", p.name)
        if m:
            names.append(m.group(1))
    seen, out = set(), []
    for n in names:
        if n not in seen:
            out.append(n); seen.add(n)
    if not out and (Path("output.pdb").exists() or Path("ligand.pdb").exists()):
        out = ["ligand"]
    return out

def _discover_jobs():
    """
    Returns a list of jobs:
      { "name": N, "orig_pdb": N.pdb, "stripped_pdbqt": N.pdbqt, "docked_pdb": output_N.pdb, "template_mol2": N.mol2 }
    with robust fallbacks for legacy single-ligand naming.
    """
    names = _read_names_from_index() or _fallback_discover_names()
    jobs = []
    single_mode = (len(names) == 1 and names[0] == "ligand")

    for nm in names:
        orig = Path(f"{nm}.pdb")
        pq   = Path(f"{nm}.pdbqt")
        mol2 = Path(f"{nm}.mol2")
        dock = Path(f"output_{nm}.pdb")
        if not dock.exists():
            bd = Path(f"best_docked_{nm}.pdb")
            dock = bd if bd.exists() else dock

        # legacy single-ligand fallbacks
        if single_mode:
            if not orig.exists() and Path("ligand.pdb").exists():     orig = Path("ligand.pdb")
            if not pq.exists()   and Path("ligand.pdbqt").exists():   pq   = Path("ligand.pdbqt")
            if not mol2.exists() and Path("ligand.mol2").exists():    mol2 = Path("ligand.mol2")
            if not dock.exists() and Path("output.pdb").exists():     dock = Path("output.pdb")

        missing = [k for k,p in {
            "orig_pdb":orig, "stripped_pdbqt":pq, "docked_pdb":dock, "template_mol2":mol2
        }.items() if not p.exists()]
        if missing:
            print(f"⚠️  Skipping {nm}: missing files {missing}")
            continue

        jobs.append({
            "name": nm, "orig_pdb": orig, "stripped_pdbqt": pq,
            "docked_pdb": dock, "template_mol2": mol2
        })

    if not jobs:
        raise SystemExit("No ligand triplets found for alignment. "
                         "Check ligands.index and Step 3 outputs (output_<name>.pdb).")
    return jobs

def _write_sdf_via_obabel(fixed_pdb: Path, out_sdf: Path):
    # Use OBabel exactly like your original pipeline
    # (this route handled phosphate O1- nicely in your prior study)
    subprocess.run(["obabel", str(fixed_pdb), "-O", str(out_sdf)])

def _write_sdf_fallback_rdkit(tag: str, mol2_path: Path, coords_xyz: np.ndarray, out_sdf: Path):
    # Very forgiving RDKit fallback: sanitize=False, then just write the coordinates to an SDF.
    from rdkit import Chem
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
    if mol is None:
        raise RuntimeError(f"[{tag}] RDKit failed to read MOL2 even with sanitize=False: {mol2_path}")
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
    print(f"   (fallback) Wrote {out_sdf} via RDKit (sanitize=False).")

def _align_one(job):
    nm = job["name"]

    docked_coords, _ = extract_coordinates(str(job["docked_pdb"]))
    full_coords, full_atoms = extract_coordinates(str(job["orig_pdb"]))
    stripped_coords, _ = extract_coordinates(str(job["stripped_pdbqt"]))

    print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms for {nm}")

    # Stabilize pairing: sort by radial distance to centroid
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned = (full_coords - c_stripped) @ R.T + c_docked  # Å

    aligned_pdb = Path(f"aligned_{nm}_fixed.pdb")
    fixed_pdb   = Path(f"fixed_{nm}.pdb")
    out_sdf     = Path(f"{nm}.sdf")

    # PDB for QC (aligned geometry)
    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")

    # Fix element column for OpenMM sanity
    fix_pdb_element_column(str(aligned_pdb), str(fixed_pdb))
    print(f"✅ Fixed PDB saved as {fixed_pdb}")

    # Write SDF via OBabel; if that fails, fall back to RDKit (sanitize=False)
    try:
        _write_sdf_via_obabel(fixed_pdb, out_sdf)
    except subprocess.CalledProcessError as e:
        print(f"⚠️  OpenBabel failed for {nm} (will try RDKit fallback).")
        try:
            _write_sdf_fallback_rdkit(nm, job["template_mol2"], aligned, out_sdf)
        except Exception as ee:
            detail = e.stderr.decode() if getattr(e, "stderr", None) else str(e)
            raise SystemExit(f"[{nm}] Could not write SDF.\n  OBabel error: {detail}\n  RDKit fallback error: {ee}")

    print(f"✅ {nm}: aligned → {aligned_pdb} | fixed → {fixed_pdb} | SDF → {out_sdf}")

def main():
    jobs = _discover_jobs()
    print(f"\n[4] Aligning {len(jobs)} ligand(s): {', '.join(j['name'] for j in jobs)}\n")
    for j in jobs:
        _align_one(j)

    # Index of SDFs for Step 5 (and a friendly summary)
    names = [j["name"] for j in jobs]
    Path("ligands_sdf.index").write_text("\n".join(names) + "\n", encoding="utf-8")
    print("\nSummary:")
    print(" • ligands_sdf.index (one name per line, SDF = <name>.sdf)")
    for n in names:
        print(f"   - {n}.sdf")

if __name__ == "__main__":
    main()
