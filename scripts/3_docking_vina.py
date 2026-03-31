# 3_docking_vina.py
# Author: Iori Mochizuki + patch
# Updated: 2025-12-02
#
# Description:
#   Run docking with AutoDock Vina for one or many ligands.
#   - Primary discovery uses 'ligands.index' written by Step 1
#   - NEW: --ligand-select allows restricting which ligands to dock
#          (by 1-based index or by name)
#   - Output per ligand: output_<name>.pdbqt, best_docked_<name>.pdb, output_<name>.pdb

import argparse, json, re, subprocess
from pathlib import Path

VINA_BIN = "./vina_1.2.7_linux_x86_64"
RECEPTOR_PDBQT = Path("receptor.pdbqt")
RECEPTOR_FOR_CENTER = Path("receptor_for_centroid.pdb")

NAME_SAFE = re.compile(r"[^A-Za-z0-9_.-]+")

def _safe(s: str) -> str:
    return NAME_SAFE.sub("_", s.strip()) or "ligand"

def _read_names_from_index():
    p = Path("ligands.index")
    if not p.exists():
        return None
    names = [ln.strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
    return [_safe(nm) for nm in names] if names else None

def _read_names_from_json():
    p = Path("ligands.json")
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        return [_safe(k) for k in d.keys()]
    except Exception:
        return None

def _discover_pdbqt_files():
    """Fallback: find *.pdbqt that look like ligands (exclude receptor & outputs)."""
    bad = {"receptor.pdbqt"}
    out = []
    for q in sorted(Path(".").glob("*.pdbqt")):
        n = q.name
        if n in bad or n.startswith("output_") or n.startswith("best_docked_"):
            continue
        out.append(q)
    return out

def _resolve_ligands():
    """
    Returns list of (name, path_to_pdbqt) in intended order.
    Primary: ligands.index → <name>.pdbqt
    Fallback: ligands.json  → <name>.pdbqt
    Last: glob *.pdbqt (excluding receptor & outputs)
    """
    # 1) ligands.index
    names = _read_names_from_index()
    if names:
        pairs = []
        missing = []
        for nm in names:
            p = Path(f"{nm}.pdbqt")
            if p.exists():
                pairs.append((nm, p))
            else:
                missing.append(nm)
        if missing:
            raise SystemExit(
                "Some ligands listed in ligands.index have no matching *.pdbqt:\n  "
                + "\n  ".join(f"- {m}.pdbqt (missing)" for m in missing)
                + "\nDid you run Step 1 in the same folder?"
            )
        return pairs

    # 2) ligands.json
    names = _read_names_from_json()
    if names:
        pairs = []
        missing = []
        for nm in names:
            p = Path(f"{nm}.pdbqt")
            if p.exists():
                pairs.append((nm, p))
            else:
                missing.append(nm)
        if pairs:
            if missing:
                print("⚠️  Warning: these entries from ligands.json were not found as *.pdbqt:")
                for m in missing:
                    print("   -", m)
            return pairs

    # 3) pure glob fallback
    files = _discover_pdbqt_files()
    if not files:
        raise SystemExit("No ligand PDBQT files found. Expected <name>.pdbqt or ligand.pdbqt after Step 1.")
    pairs = [(_safe(f.stem), f) for f in files]
    return pairs

def _filter_ligands_by_selector(all_pairs, selector: str | None):
    """
    all_pairs: list of (name, Path)
    selector: comma-separated list of 1-based indices or names (stems)
      - "1,3" → first and third in all_pairs
      - "lucidenic_acid_A_reishi_,procyanidin_C1_GSE_" → by name
    If selector is None, returns all_pairs unchanged.
    """
    if not selector:
        return all_pairs

    tokens = [t.strip() for t in selector.split(",") if t.strip()]
    if not tokens:
        return all_pairs

    # 1-based index map
    by_idx = {str(i): all_pairs[i-1] for i in range(1, len(all_pairs)+1)}
    # name map
    by_name = {nm: (nm, p) for nm, p in all_pairs}

    selected = []
    for tok in tokens:
        if tok in by_idx:
            selected.append(by_idx[tok])
        elif tok in by_name:
            selected.append(by_name[tok])
        else:
            print(f"⚠️  --ligand-select token '{tok}' not matched; ignoring.")

    # Keep order consistent with all_pairs
    seen = set()
    ordered = []
    for nm, p in all_pairs:
        if (nm, p) in selected and (nm, p) not in seen:
            ordered.append((nm, p))
            seen.add((nm, p))

    if not ordered:
        print("⚠️  --ligand-select did not match any ligands; using all.")
        return all_pairs
    return ordered

def _compute_center(args):
    # Priority 1: manual center
    if args.center_x is not None and args.center_y is not None and args.center_z is not None:
        x = round(float(args.center_x), 3)
        y = round(float(args.center_y), 3)
        z = round(float(args.center_z), 3)
        print(f"📍 Docking box center (manual): x={x} y={y} z={z}")
        return x, y, z

    # Priority 2: RNA residue range centroid
    def _parse_range(s):
        for dash in ['–', '—', '−', '‒', '―']:
            s = s.replace(dash, '-')
        a, b = s.split('-')
        return int(a), int(b)

    txt = RECEPTOR_FOR_CENTER.read_text()  # raises if missing (better fail fast)

    if args.rna and args.res_range:
        start, end = _parse_range(args.res_range)
        sx = sy = sz = 0.0
        n = 0
        for ln in txt.splitlines():
            if ln.startswith("ATOM"):
                resid = int(ln[22:26])
                if start <= resid <= end:
                    sx += float(ln[30:38])
                    sy += float(ln[38:46])
                    sz += float(ln[46:54])
                    n += 1
        if n > 0:
            x, y, z = round(sx/n, 3), round(sy/n, 3), round(sz/n, 3)
            print(f"📍 RNA centroid {start}-{end}: x={x} y={y} z={z}")
            return x, y, z
        print("⚠️  No atoms in the specified RNA range; falling back to HETATM/ATOM centroid.")

    # Priority 3a: specific HETATM names (e.g., 3E4)
    if args.het_resnames:
        wanted = {t.strip().upper() for t in args.het_resnames.split(",") if t.strip()}
        sx = sy = sz = 0.0
        n = 0
        for ln in txt.splitlines():
            if ln.startswith("HETATM") and ln[17:20].strip().upper() in wanted:
                sx += float(ln[30:38])
                sy += float(ln[38:46])
                sz += float(ln[46:54])
                n += 1
        if n > 0:
            x, y, z = round(sx/n, 3), round(sy/n, 3), round(sz/n, 3)
            print(f"📍 HETATM({','.join(sorted(wanted))}) centroid: x={x} y={y} z={z}")
            return x, y, z
        print(f"⚠️  No HETATM with names in {wanted}; falling back.")

    # Priority 3b: HETATM centroid (if requested)
    if args.use_residue_centroid:
        resnames = set()
        for ln in txt.splitlines():
            if ln.startswith("HETATM"):
                resnames.add(ln[17:20].strip())
        sx = sy = sz = 0.0
        n = 0
        for ln in txt.splitlines():
            if ln.startswith("HETATM") and (ln[17:20].strip() in resnames):
                sx += float(ln[30:38])
                sy += float(ln[38:46])
                sz += float(ln[46:54])
                n += 1
        if n > 0:
            x, y, z = round(sx/n, 3), round(sy/n, 3), round(sz/n, 3)
            print(f"📍 HETATM centroid: x={x} y={y} z={z}")
            return x, y, z
        print("⚠️  No HETATM residues; falling back to ATOM centroid.")

    # Priority 4: all ATOM centroid
    sx = sy = sz = 0.0
    n = 0
    for ln in txt.splitlines():
        if ln.startswith("ATOM"):
            sx += float(ln[30:38])
            sy += float(ln[38:46])
            sz += float(ln[46:54])
            n += 1
    x, y, z = round(sx/n, 3), round(sy/n, 3), round(sz/n, 3)
    print(f"📍 ATOM centroid: x={x} y={y} z={z}")
    return x, y, z

def _extract_model1_pdb(in_pdbqt: Path, out_pdb: Path):
    """Write MODEL 1’s (HET)ATOM/TER records to out_pdb."""
    inside = False
    with in_pdbqt.open("r") as fin, out_pdb.open("w") as fout:
        for line in fin:
            if line.startswith("MODEL 2"):
                break
            if line.startswith("MODEL 1"):
                inside = True
                continue
            if not inside:
                continue
            if line.startswith(("ATOM", "HETATM", "TER")):
                fout.write(line)

def main():
    ap = argparse.ArgumentParser(description="Dock one or many ligands with AutoDock Vina.")
    ap.add_argument("--use-residue-centroid", action="store_true",
                    help="Use centroid of all HETATM residues (unless overridden by manual center or het_resnames).")
    ap.add_argument("--rna", action="store_true", help="Enable RNA centroid mode.")
    ap.add_argument("--res-range", type=str, default=None,
                    help="Residue range for RNA centroid, e.g. '10-50'.")
    ap.add_argument("--center_x", type=float, default=None)
    ap.add_argument("--center_y", type=float, default=None)
    ap.add_argument("--center_z", type=float, default=None)
    ap.add_argument("--vina_seed", type=int, default=12345)
    ap.add_argument("--het-resnames", type=str, default=None,
                    help="Comma-separated HETATM 3-letter names to centroid (e.g., '3E4,EPE').")
    ap.add_argument("--ligand-select", type=str, default=None,
                    help="Subset of ligands to dock, by 1-based index or name (comma separated).")

    args = ap.parse_args()

    if not RECEPTOR_PDBQT.exists():
        raise SystemExit("Missing receptor.pdbqt in current directory.")
    if not RECEPTOR_FOR_CENTER.exists():
        raise SystemExit("Missing receptor_for_centroid.pdb in current directory.")

    ligands = _resolve_ligands()
    print(f"Found {len(ligands)} ligand(s): " + ", ".join(nm for nm, _ in ligands))

    # NEW: apply selector, just like Step 5
    ligands = _filter_ligands_by_selector(ligands, args.ligand_select)
    print("Docking the following ligand(s): " + ", ".join(nm for nm, _ in ligands))

    cx, cy, cz = _compute_center(args)

    legacy_single = (len(ligands) == 1 and ligands[0][1].name == "ligand.pdbqt")

    for i, (name, ligp) in enumerate(ligands, start=1):
        tag = _safe(name)
        out_pdbqt = Path(f"output_{tag}.pdbqt")
        out_best  = Path(f"best_docked_{tag}.pdb")
        out_clean = Path(f"output_{tag}.pdb")

        cmd = (
            f"{VINA_BIN} "
            f"--receptor {RECEPTOR_PDBQT} "
            f"--ligand {ligp} "
            f"--out {out_pdbqt} "
            f"--center_x {cx} --center_y {cy} --center_z {cz} "
            f"--size_x 30 --size_y 30 --size_z 30 "
            f"--seed {int(args.vina_seed)} --exhaustiveness 20"
        )
        print(f"\n▶️  [{i}/{len(ligands)}] Vina docking for {ligp.name}")
        subprocess.run(cmd, shell=True, check=True)
        _extract_model1_pdb(out_pdbqt, out_best)
        out_clean.write_text(out_best.read_text())
        print(f"✅ {tag}: best pose → {out_best}  | cleaned → {out_clean}")

        if legacy_single:
            Path("output.pdbqt").write_bytes(out_pdbqt.read_bytes())
            Path("best_docked_ligand.pdb").write_bytes(out_best.read_bytes())
            Path("output.pdb").write_bytes(out_clean.read_bytes())
            print("↪️  Also wrote legacy output.pdbqt / best_docked_ligand.pdb / output.pdb")

if __name__ == "__main__":
    main()
