# 4_align_ligand.py
# Author: Iori Mochizuki
# Updated: 2025-10-30
# Description: Align many ligands to their docked coordinates; write per‑ligand PDB+SDF.
# Policy: SDF is written with Open Babel from the fixed aligned PDB (no RDKit fallbacks).

from pathlib import Path
import re
import subprocess
import sys
import numpy as np

from utils import extract_coordinates, kabsch, fix_pdb_element_column

# Step 1 writes this; we consume it here.
INDEX_FILE = Path("ligands.index")

def _read_names_from_index() -> list[str]:
    """Read unique, non-comment ligand names from ligands.index."""
    if INDEX_FILE.exists():
        names = []
        seen = set()
        for ln in INDEX_FILE.read_text(encoding="utf-8").splitlines():
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            if s not in seen:
                names.append(s)
                seen.add(s)
        return names
    return []

def _fallback_discover_names() -> list[str]:
    """
    If no index is present, derive names from Step 3 products:
      output_<name>.pdb → "<name>"
    Keep order by filename; deduplicate.
    Legacy single-ligand fallback → ["ligand"] if needed.
    """
    out = []
    seen = set()
    for p in sorted(Path(".").glob("output_*.pdb")):
        m = re.match(r"output_(.+)\.pdb$", p.name)
        if m:
            nm = m.group(1)
            if nm not in seen:
                out.append(nm)
                seen.add(nm)
    if not out and (Path("output.pdb").exists() or Path("ligand.pdb").exists()):
        out = ["ligand"]
    return out

def _discover_jobs():
    """
    Build per‑ligand job dicts:
      {
        "name": <nm>,
        "orig_pdb": <nm>.pdb (or ligand.pdb),
        "stripped_pdbqt": <nm>.pdbqt (or ligand.pdbqt),
        "docked_pdb": output_<nm>.pdb (or best_docked_<nm>.pdb; legacy: output.pdb)
      }
    """
    names = _read_names_from_index() or _fallback_discover_names()
    if not names:
        raise SystemExit("No ligands to align. Provide ligands.index or ensure output_<name>.pdb exists from Step 3.")

    jobs = []
    single_mode = (len(names) == 1 and names[0] == "ligand")
    for nm in names:
        orig = Path(f"{nm}.pdb")
        pq   = Path(f"{nm}.pdbqt")
        dock = Path(f"output_{nm}.pdb")
        if not dock.exists():
            alt = Path(f"best_docked_{nm}.pdb")
            if alt.exists():
                dock = alt

        # Legacy single‑ligand names
        if single_mode:
            if not orig.exists() and Path("ligand.pdb").exists():     orig = Path("ligand.pdb")
            if not pq.exists()   and Path("ligand.pdbqt").exists():   pq   = Path("ligand.pdbqt")
            if not dock.exists() and Path("output.pdb").exists():     dock = Path("output.pdb")

        missing = [k for k, p in {
            "orig_pdb": orig, "stripped_pdbqt": pq, "docked_pdb": dock
        }.items() if not p.exists()]
        if missing:
            print(f"⚠️  Skipping {nm}: missing files {missing}")
            continue

        jobs.append({
            "name": nm,
            "orig_pdb": orig,
            "stripped_pdbqt": pq,
            "docked_pdb": dock
        })

    if not jobs:
        raise SystemExit("No complete ligand triplets found (orig PDB, stripped PDBQT, docked PDB).")
    return jobs

def _obabel_pdb_to_sdf(fixed_pdb: Path, out_sdf: Path):
    """Use Open Babel exactly as in the original single‑ligand pipeline."""
    try:
        subprocess.run(["obabel", str(fixed_pdb), "-O", str(out_sdf)],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise SystemExit("Open Babel ('obabel') not found in PATH. Install it in your conda env.")
    except subprocess.CalledProcessError as e:
        err = e.stderr.decode(errors="ignore") if e.stderr else str(e)
        raise SystemExit(f"Open Babel failed converting {fixed_pdb.name} → {out_sdf.name}:\n{err}")

def _align_one(job):
    nm = job["name"]

    docked_coords, _ = extract_coordinates(str(job["docked_pdb"]))
    full_coords, full_atoms = extract_coordinates(str(job["orig_pdb"]))
    stripped_coords, _ = extract_coordinates(str(job["stripped_pdbqt"]))

    print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms for {nm}")

    # Stabilize pairing: sort by radial distance to each centroid
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    # Kabsch alignment (utils.kabsch)
    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned = (full_coords - c_stripped) @ R.T + c_docked  # Å

    aligned_pdb = Path(f"aligned_{nm}_fixed.pdb")
    fixed_pdb   = Path(f"fixed_{nm}.pdb")
    out_sdf     = Path(f"{nm}.sdf")

    # Write aligned PDB (raw coordinates applied onto original atom order/records)
    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")

    # Fix element column for OpenMM sanity (same helper as single‑ligand flow)
    fix_pdb_element_column(str(aligned_pdb), str(fixed_pdb))
    print(f"✅ Fixed PDB saved as {fixed_pdb}")

    # Convert to SDF via Open Babel (no RDKit fallback by design)
    _obabel_pdb_to_sdf(fixed_pdb, out_sdf)

    print(f"✅ {nm}: aligned → {aligned_pdb} | fixed → {fixed_pdb} | SDF → {out_sdf}")

def main():
    jobs = _discover_jobs()
    print(f"\n[4] Aligning {len(jobs)} ligand(s): {', '.join(j['name'] for j in jobs)}\n")

    for j in jobs:
        _align_one(j)

    # Write the SDF index for Step 5 (deterministic order = names order here)
    names = [j["name"] for j in jobs]
    Path("ligands_sdf.index").write_text("\n".join(names) + "\n", encoding="utf-8")

    print("\nSummary:")
    print(" • ligands_sdf.index (one name per line, SDF = <name>.sdf)")
    for n in names:
        print(f"   - {n}.sdf")

if __name__ == "__main__":
    main()
