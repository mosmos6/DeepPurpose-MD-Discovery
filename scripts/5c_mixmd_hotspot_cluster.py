# -*- coding: utf-8 -*-
"""
Greedy hotspot clustering for MixMD probe centroids (dependency‑free).

Input  : CSV from ProbeHotspotReporter (columns: step,resname,resid,x_nm,y_nm,z_nm)
Method : For each probe type separately:
           1) time‑average (x,y,z) per unique probe instance (resname,resid)
           2) greedy clustering with merge radius R (nm) on the averaged points
           3) rank clusters by occupancy (unique instances captured)
Outputs:
  - build/hotspots_<stem>_r{R}nm.csv     (ranked cluster table)
  - build/hotspots_<stem>_r{R}nm.pdb     (one pseudo‑atom per cluster center)
  - build/hotspots_<stem>_r{R}nm.vina.txt  (optional docking box hints)

Usage example:
  python scripts/5c_mixmd_hotspot_cluster.py \
      --input mixmd_hotspots_no_ligand.csv \
      --radius-nm 0.35 \
      --min-members 3 \
      --top-k 3
"""

import argparse, csv, math
from pathlib import Path

# --- Resname handling (accept both 3‑letter canonical and 4‑letter aliases) ---
CANONICAL = {
    "IPA": "IPA", "ACN": "ACN", "IMD": "IMD",
    "ACM": "ACM", "PHO": "PHO", "HAC": "HAC",
    # aliases accepted in the CSV (will be canonicalized below)
    "ACEA": "ACM", "PHOL": "PHO", "ACOH": "HAC"
}
# Pick a consistent order for output
PROBE_ORDER = ["IPA", "ACN", "IMD", "ACM", "PHO", "HAC"]

def canon_name(s: str) -> str:
    s = (s or "").strip().upper()
    return CANONICAL.get(s, s)

# --- Small geometry helpers (no NumPy) ---
def d2(a, b) -> float:
    """Squared distance in nm^2 between 3‑tuples a and b."""
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return dx*dx + dy*dy + dz*dz

def mean_of(points):
    """Plain average of 3‑tuples."""
    n = len(points)
    sx = sy = sz = 0.0
    for (x,y,z) in points:
        sx += x; sy += y; sz += z
    return (sx/n, sy/n, sz/n)

# --- Core: greedy clustering on averaged instance positions ---
def greedy_clusters(points, radius_nm: float, refine_iters: int = 1):
    """
    points: list[(x,y,z)] in nm (already time‑averaged per instance)
    radius_nm: merge radius (nm)
    refine_iters: optional small center‑recompute passes
    Returns list of dicts: {"center":(x,y,z),"members":[indices]}
    """
    R2 = radius_nm * radius_nm
    n = len(points)
    unassigned = [True]*n

    # precompute neighbor counts to start with densest point
    neighbor_counts = [0]*n
    for i in range(n):
        pi = points[i]
        c = 0
        for j in range(n):
            if d2(pi, points[j]) <= R2:
                c += 1
        neighbor_counts[i] = c

    clusters = []
    remaining = sum(unassigned)
    while remaining > 0:
        # pick densest remaining
        seed = -1
        best = -1
        for i in range(n):
            if unassigned[i] and neighbor_counts[i] > best:
                best = neighbor_counts[i]; seed = i
        if seed < 0:
            break

        # initial membership from seed
        center = points[seed]
        members = []
        for i in range(n):
            if unassigned[i] and d2(points[i], center) <= R2:
                members.append(i)

        # refine center/membership a few times
        for _ in range(refine_iters):
            if not members:
                break
            center = mean_of([points[i] for i in members])
            new_members = []
            for i in range(n):
                if unassigned[i] and d2(points[i], center) <= R2:
                    new_members.append(i)
            if set(new_members) == set(members):
                break
            members = new_members

            # safety: if refinement drifted, keep it bounded
            if not members:
                members = [seed]; center = points[seed]; break

        # finalize this cluster
        clusters.append({"center": center, "members": members})
        for i in members:
            unassigned[i] = False
        remaining = sum(unassigned)

    # sort by descending size
    clusters.sort(key=lambda c: len(c["members"]), reverse=True)
    return clusters

# --- I/O helpers ---
def read_hotspot_csv(path: Path):
    """Return rows = list[{step:int,resname:str,resid:int,x_nm:float,y_nm:float,z_nm:float}]"""
    rows = []
    with open(path, newline="") as f:
        # skip commented diagnostics/header lines (starting with '#')
        first = f.readline()
        if not first.startswith("#"):
            # rewind if first line already header
            f.seek(0)
        r = csv.DictReader(f)
        for rec in r:
            try:
                rows.append({
                    "step":   int(rec["step"]),
                    "resname": canon_name(rec["resname"]),
                    "resid":  int(rec["resid"]),
                    "x_nm":   float(rec["x_nm"]),
                    "y_nm":   float(rec["y_nm"]),
                    "z_nm":   float(rec["z_nm"]),
                })
            except Exception:
                # tolerate incomplete rows
                continue
    return rows

def time_average_by_instance(rows):
    """
    Group by (resname,resid) and return per‑probe:
      dict[probe] = list of {"avg":(x,y,z), "count":n, "key":(resname,resid)}
    """
    acc = {}  # (resname,resid) -> [n, sx, sy, sz]
    for r in rows:
        k = (r["resname"], r["resid"])
        e = acc.get(k)
        if e is None:
            acc[k] = [0, 0.0, 0.0, 0.0]
            e = acc[k]
        e[0] += 1; e[1] += r["x_nm"]; e[2] += r["y_nm"]; e[3] += r["z_nm"]

    out = {}
    for (resname, resid), (n, sx, sy, sz) in acc.items():
        if n <= 0:
            continue
        avg = (sx/n, sy/n, sz/n)
        out.setdefault(resname, []).append({"avg": avg, "count": n, "key": (resname, resid)})
    return out

def ensure_build_dir():
    b = Path("build"); b.mkdir(exist_ok=True, parents=True); return b

def write_clusters_csv(out_csv: Path, clusters_by_probe, radius_nm, top_k=None, min_members=1):
    """
    clusters_by_probe: dict[probe] = list of {"center":(x,y,z),"members":[idx], "instances":[(resname,resid),...]}
    Writes one table with Å & nm coordinates and occupancy.
    """
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "probe","cluster_rank","count",
            "x_nm","y_nm","z_nm",
            "x_A","y_A","z_A"
        ])
        for probe in PROBE_ORDER:
            clusters = clusters_by_probe.get(probe, [])
            # occupancy max for normalization (for PDB occupancy later)
            # but here we just output centers
            kept = 0
            for rank, c in enumerate(clusters, start=1):
                cnt = len(c["members"])
                if cnt < min_members:
                    continue
                x,y,z = c["center"]
                w.writerow([probe, rank, cnt, f"{x:.5f}", f"{y:.5f}", f"{z:.5f}",
                            f"{x*10:.3f}", f"{y*10:.3f}", f"{z*10:.3f}"])
                kept += 1
                if top_k and kept >= top_k:
                    break

def write_clusters_pdb(out_pdb: Path, clusters_by_probe, top_k=None, min_members=1):
    """
    Write one pseudo‑atom per cluster:
      - resname = probe
      - occupancy = count / max_count_for_that_probe
      - B‑factor  = count (raw)
    """
    lines = []
    serial = 1
    # fixed chain map for readability
    chain_for = {"IPA":"I","ACN":"N","IMD":"M","ACM":"A","PHO":"P","HAC":"H"}
    for probe in PROBE_ORDER:
        clusters = clusters_by_probe.get(probe, [])
        if not clusters:
            continue
        maxcnt = max(len(c["members"]) for c in clusters) if clusters else 1

        kept = 0
        for rank, c in enumerate(clusters, start=1):
            cnt = len(c["members"])
            if cnt < min_members:
                continue
            x,y,z = c["center"]
            xA,yA,zA = x*10.0, y*10.0, z*10.0
            occ = min(1.0, cnt / float(maxcnt))
            bfac = float(cnt)

            atom_name = "HSP "  # generic pseudo atom name
            resname = probe[:3].rjust(3)
            chain = chain_for.get(probe, "X")
            resid = rank  # rank within that probe

            lines.append(
                f"HETATM{serial:5d} {atom_name:4s}{resname:>4s}{chain}{resid:4d}    "
                f"{xA:8.3f}{yA:8.3f}{zA:8.3f}{occ:6.2f}{bfac:6.2f}          {'X':>2s}"
            )
            serial += 1
            kept += 1
            if top_k and kept >= top_k:
                break

    lines.append("END\n")
    out_pdb.write_text("\n".join(lines), encoding="utf-8")

def write_vina_boxes(out_txt: Path, clusters_by_probe, size_A=24.0, top_k=3, min_members=1):
    """
    Small helper: one line per cluster usable as Vina box snippet.
    center_(x|y|z) in Å; size_(x|y|z) = size_A (cubic).
    """
    with open(out_txt, "w") as f:
        f.write("# AutoDock Vina box suggestions (Å)\n")
        f.write("# columns: probe  rank  count  center_x  center_y  center_z  size\n")
        for probe in PROBE_ORDER:
            clusters = clusters_by_probe.get(probe, [])
            kept = 0
            for rank, c in enumerate(clusters, start=1):
                cnt = len(c["members"])
                if cnt < min_members:
                    continue
                xA,yA,zA = (c["center"][0]*10.0, c["center"][1]*10.0, c["center"][2]*10.0)
                f.write(f"{probe:>4s}  {rank:2d}  {cnt:4d}  {xA:8.3f} {yA:8.3f} {zA:8.3f}  {size_A:.1f}\n")
                kept += 1
                if top_k and kept >= top_k:
                    break

def main():
    ap = argparse.ArgumentParser(description="Greedy hotspot clustering for MixMD probe centroids.")
    ap.add_argument("--input", required=True, help="Hotspot CSV from 5_md_simulation.py reporter.")
    ap.add_argument("--radius-nm", type=float, default=0.35, help="Merge radius (nm). ~0.30–0.40 nm works well.")
    ap.add_argument("--refine-iters", type=int, default=1, help="Center refinement passes per cluster.")
    ap.add_argument("--min-members", type=int, default=2, help="Keep clusters with at least this many distinct instances.")
    ap.add_argument("--top-k", type=int, default=3, help="Write only top‑K clusters per probe to PDB/boxes (0 = all).")
    ap.add_argument("--vina-size-A", type=float, default=24.0, help="Suggested cubic box edge for Vina (Å).")
    args = ap.parse_args()

    inp = Path(args.input)
    rows = read_hotspot_csv(inp)
    if not rows:
        raise SystemExit(f"No rows parsed from {inp}")

    # 1) time‑average per (resname,resid)
    by_probe = time_average_by_instance(rows)

    # 2) greedy clustering per probe
    clusters_by_probe = {}
    for probe in PROBE_ORDER:
        instances = by_probe.get(probe, [])
        if not instances:
            continue
        pts = [rec["avg"] for rec in instances]
        cls = greedy_clusters(pts, radius_nm=args.radius_nm, refine_iters=args.refine_iters)
        clusters_by_probe[probe] = cls

    # 3) outputs
    build = ensure_build_dir()
    stem = inp.stem  # e.g., mixmd_hotspots_no_ligand
    tag = f"{stem}_r{args.radius_nm:.2f}nm".replace(".", "p")
    out_csv = build / f"hotspots_{tag}.csv"
    out_pdb = build / f"hotspots_{tag}.pdb"
    out_txt = build / f"hotspots_{tag}.vina.txt"

    write_clusters_csv(out_csv, clusters_by_probe, args.radius_nm, top_k=args.top_k or None, min_members=args.min_members)
    write_clusters_pdb(out_pdb, clusters_by_probe, top_k=args.top_k or None, min_members=args.min_members)
    write_vina_boxes(out_txt, clusters_by_probe, size_A=float(args.vina_size_A), top_k=args.top_k or None, min_members=args.min_members)

    # Small console summary
    print(f"[6] Wrote clusters → {out_csv}")
    print(f"[6] Centers PDB    → {out_pdb}")
    print(f"[6] Vina boxes     → {out_txt}")

if __name__ == "__main__":
    main()
