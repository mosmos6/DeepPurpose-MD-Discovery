# 5c_mixmd_hotspot_cluster.py
# Greedy clustering of MixMD probe hotspots (dependency-free).
# Reads centroids CSV from 5_md_simulation.py's ProbeHotspotReporter.

import argparse, csv, math, sys
from collections import defaultdict, Counter
from pathlib import Path

# ---------- canonical names & aliases ----------
ALIASES_TO_CANON = {
    # canonical (3-letter)
    "IPA": "IPA", "ACN": "ACN", "IMD": "IMD", "ACM": "ACM", "PHO": "PHO", "HAC": "HAC",
    # legacy aliases used in older notes/CLIs
    "ACEA": "ACM", "PHOL": "PHO", "ACOH": "HAC",
}
CANON_SET = {"IPA","ACN","IMD","ACM","PHO","HAC"}

# nice, distinct chain IDs per type (single char; viewer-friendly)
CHAIN_ID = {"IPA":"I","ACN":"A","IMD":"M","ACM":"C","PHO":"P","HAC":"H"}

# ---------- I/O helpers ----------
def _unify_resname(name: str) -> str | None:
    if not name:
        return None
    key = name.strip().upper()
    return ALIASES_TO_CANON.get(key, None)

def _iter_hotspot_rows(csv_path: Path):
    """
    Yield (step:int, resname:canonical, resid:int, x:float, y:float, z:float)
    Skips comments (lines starting with '#').
    Accepts headers with x_nm/y_nm/z_nm (preferred) or x/y/z.
    """
    with open(csv_path, "r", newline="") as f:
        rows = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    if not rows:
        return

    rdr = csv.DictReader(rows)

    def _get(row, *keys, default=None):
        for k in keys:
            if k in row:
                return row[k]
        # lowercase fallback
        low = {k.lower(): v for k, v in row.items()}
        for k in keys:
            if k.lower() in low:
                return low[k.lower()]
        return default

    for row in rdr:
        rn_raw = _get(row, "resname")
        rn = _unify_resname(rn_raw)
        if rn is None:
            continue
        try:
            step  = int(float(_get(row, "step", default="0")))
            resid = int(float(_get(row, "resid", default="0")))
            x = float(_get(row, "x_nm", "x"))
            y = float(_get(row, "y_nm", "y"))
            z = float(_get(row, "z_nm", "z"))
        except Exception:
            continue
        yield step, rn, resid, x, y, z

# ---------- core logic ----------
def _time_average_by_instance(rows, min_frames: int):
    """
    rows: iterable of (step, rn, resid, x,y,z)
    Return dict[(rn,resid)] -> (rn, resid, n_frames, x̄,ȳ,ż)
    """
    acc = {}
    for _, rn, resid, x, y, z in rows:
        key = (rn, resid)
        if key not in acc:
            acc[key] = [0.0, 0.0, 0.0, 0]  # sx, sy, sz, n
        acc[key][0] += x
        acc[key][1] += y
        acc[key][2] += z
        acc[key][3] += 1

    out = {}
    for (rn, resid), (sx, sy, sz, n) in acc.items():
        if n >= min_frames:
            out[(rn, resid)] = (rn, resid, n, sx/n, sy/n, sz/n)
    return out

def _dist(a, b):
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def _greedy_cluster(points, radius_nm):
    """
    points: list of (rn, resid, n_frames, x, y, z) FOR A SINGLE rn (type)
    Returns list of clusters:
      { 'type': rn, 'centroid': (x,y,z), 'members': [(resid,n_frames), ...] }
    Ranking: by number of distinct members (occupancy).
    """
    clusters = []
    for rn, resid, nfr, x, y, z in points:
        p = (x, y, z)
        assigned = False
        for cl in clusters:
            if _dist(p, cl['centroid']) <= radius_nm:
                m = len(cl['members'])
                cx, cy, cz = cl['centroid']
                cl['centroid'] = ((cx*m + x)/(m+1), (cy*m + y)/(m+1), (cz*m + z)/(m+1))
                cl['members'].append((resid, nfr))
                assigned = True
                break
        if not assigned:
            clusters.append({'type': rn, 'centroid': p, 'members': [(resid, nfr)]})
    clusters.sort(key=lambda c: len(c['members']), reverse=True)
    return clusters

# ---------- PDB writer ----------
def _write_pdb(base_path: Path, clusters_all, append_to_pdb: Path | None = None):
    """
    Write hotspots as pseudo-atoms (HETATM) at cluster centroids.
    Coordinates are in Å (input is nm → ×10).
    If append_to_pdb is provided, copy that file first and then append hotspots.
    """
    pdb_path = base_path.with_suffix(".pdb")
    nm_to_A = 10.0

    # If appending to an existing PDB, copy content first
    lines = []
    if append_to_pdb:
        src = Path(append_to_pdb)
        if not src.exists():
            raise SystemExit(f"--append-to-pdb file not found: {src}")
        lines = src.read_text().splitlines()

        # Strip possible trailing END/ENDMDL to safely append
        while lines and lines[-1].strip() in {"END", "ENDMDL"}:
            lines.pop()

    # Prepare new HETATM records
    het_lines = []
    serial = 1
    last_chain = None
    for idx, (rn, local_rank, cl) in enumerate(clusters_all, start=1):
        cx, cy, cz = cl['centroid']
        xA, yA, zA = cx*nm_to_A, cy*nm_to_A, cz*nm_to_A
        count = len(cl['members'])
        resname = rn                 # IPA/ACN/… (3-letter)
        chain  = CHAIN_ID.get(rn, "X")
        resseq = int(local_rank)     # 1..K within that type

        # Insert TER when chain changes (helps viewers)
        if last_chain is not None and chain != last_chain:
            het_lines.append("TER".ljust(80))
        last_chain = chain

        # HETATM format (PDB v3). B-factor stores occupancy count for quick visualization.
        # Columns:  1-6 'HETATM', 7-11 serial, 13-16 atom, 17 alt, 18-20 resName,
        #           22 chain, 23-26 resSeq, 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 temp, 77-78 element
        het = (
            f"HETATM{serial:5d}  CEN {resname:>3s} {chain:1s}"
            f"{resseq:4d}    "
            f"{xA:8.3f}{yA:8.3f}{zA:8.3f}"
            f"{1.00:6.2f}{min(99.99,float(count)):6.2f}"
            f"          C "
        )
        het_lines.append(het.ljust(80))
        serial += 1

    het_lines.append("TER".ljust(80))
    het_lines.append("END".ljust(80))

    # Write final
    with open(pdb_path, "w") as f:
        if lines:
            for ln in lines:
                f.write(ln.rstrip() + "\n")
        # header remarks for the hotspot block
        f.write(f"REMARK 500 HOTSPOT CENTROIDS (per-probe clusters)\n")
        f.write(f"REMARK 500 B-FACTOR FIELD = number of probe instances in cluster\n")
        for ln in het_lines:
            f.write(ln + "\n")

    print(f"✅ Wrote: {pdb_path}")

# ---------- CLI & pipeline ----------
def main():
    ap = argparse.ArgumentParser(
        description="Cluster MixMD probe hotspots (greedy, dependency-free)."
    )
    ap.add_argument("--input", required=True,
                    help="ProbeHotspotReporter CSV, e.g., mixmd_hotspots_no_ligand.csv")
    ap.add_argument("--radius-nm", type=float, default=0.40,
                    help="Merge radius in nm (default 0.40)")
    ap.add_argument("--min-members", type=int, default=3,
                    help="Min frames per probe instance to keep (default 3)")
    ap.add_argument("--top-k", type=int, default=3,
                    help="Top K clusters per probe type to export (default 3)")
    ap.add_argument("--types", type=str, default=None,
                    help="Comma-separated probe types to include (accepts aliases). Default: all present.")
    ap.add_argument("--list-names", action="store_true",
                    help="Print residue names found in CSV and exit.")
    ap.add_argument("--append-to-pdb", type=str, default=None,
                    help="If provided, append hotspot pseudo-atoms to this PDB (useful for quick overlay).")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"Input CSV not found: {in_path}")

    # Load once for diagnostics
    all_rows = list(_iter_hotspot_rows(in_path))

    if args.list_names:
        cnt = Counter(rn for _, rn, _, _, _, _ in all_rows)
        if not cnt:
            print("No recognizable probe residues found. Names seen in file may not match expected set.")
        else:
            print("Residue names found (canonicalized):")
            for rn, n in cnt.most_common():
                print(f"  {rn}: {n} rows")
        return

    if not all_rows:
        raise SystemExit(
            f"No rows parsed from {in_path}.\n"
            f"Tip: run with --list-names to see which residue names are present."
        )

    # Optional type filter
    type_filter = None
    if args.types:
        user = {t.strip().upper() for t in args.types.split(",") if t.strip()}
        type_filter = {ALIASES_TO_CANON.get(t, t) for t in user} & set(ALIASES_TO_CANON.values())
        if not type_filter:
            print(f"Warning: none of the specified types matched known set {sorted(CANON_SET)}; using all types.")
            type_filter = None

    # Time-average per probe instance
    per_inst = _time_average_by_instance(all_rows, min_frames=max(1, args.min_members))

    # Split points by type
    points_by_type = defaultdict(list)
    for (rn, resid), (rn_, resid_, nfr, x, y, z) in per_inst.items():
        if type_filter and rn_ not in type_filter:
            continue
        points_by_type[rn_].append((rn_, resid_, nfr, x, y, z))

    if not any(points_by_type.values()):
        cnt = Counter(rn for _, rn, _, _, _, _ in all_rows)
        raise SystemExit(
            "No rows survived filtering.\n"
            f"  • Present names (counts): {dict(cnt)}\n"
            f"  • Filter types: {sorted(type_filter) if type_filter else 'None (using all)'}\n"
            "  • Check that your CSV uses IPA/ACN/IMD/ACM/PHO/HAC (aliases ACEA/PHOL/ACOH also accepted).\n"
            "  • Or try a smaller --min-members."
        )

    # Cluster each type
    clusters_all = []
    for rn, pts in points_by_type.items():
        if not pts:
            continue
        clusters = _greedy_cluster(pts, radius_nm=float(args.radius_nm))
        # keep top-k for this rn
        clusters_all.extend([(rn, i+1, c) for i, c in enumerate(clusters[:max(1, args.top_k)])])

    if not clusters_all:
        raise SystemExit("No clusters were formed. Try a larger --radius-nm or smaller --min-members.")

    # Write outputs
    base = f"hotspots_{in_path.stem}_r{float(args.radius_nm):.2f}nm"
    csv_out = Path(f"{base}.csv")
    vina_out = Path(f"{base}.vina.txt")

    with open(csv_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank","type","count","x_nm","y_nm","z_nm","member_resids","member_frames_sum"])
        rank_global = 1
        for rn, _, cl in clusters_all:
            mem = cl['members']
            count = len(mem)
            cx, cy, cz = cl['centroid']
            member_resids = ";".join(str(rid) for rid, _ in mem)
            frames_sum = sum(nfr for _, nfr in mem)
            w.writerow([rank_global, rn, count,
                        f"{cx:.3f}", f"{cy:.3f}", f"{cz:.3f}",
                        member_resids, frames_sum])
            rank_global += 1

    with open(vina_out, "w") as f:
        f.write("# AutoDock Vina boxes (center in nm; convert to Å by ×10 if needed)\n")
        for idx, (rn, _, cl) in enumerate(clusters_all, start=1):
            cx, cy, cz = cl['centroid']
            f.write(f"name=HOTSPOT_{rn}_{idx}\n")
            f.write(f"center_x_nm={cx:.4f}\ncenter_y_nm={cy:.4f}\ncenter_z_nm={cz:.4f}\n")
            f.write("size_x_A=30\nsize_y_A=30\nsize_z_A=30\n\n")

    print(f"✅ Wrote: {csv_out}")
    print(f"✅ Wrote: {vina_out}")

    # --- NEW: PDB with pseudo-atoms at cluster centroids ---
    _write_pdb(Path(base), clusters_all,
               append_to_pdb=(Path(args.append_to_pdb) if args.append_to_pdb else None))

if __name__ == "__main__":
    main()
