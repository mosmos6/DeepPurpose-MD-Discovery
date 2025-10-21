# 5c_mixmd_hotspot_cluster.py
# Greedy clustering of MixMD probe centroids with optional Kabsch alignment
# Author: Iori Mochizuki (pipeline)
# Updated: 2025-10-17

import argparse, csv, os, math
from pathlib import Path
from collections import defaultdict

# numpy is used only here (post-processing, no OpenMM units); safe & fast
import numpy as np

# -----------------------------
# I/O helpers
# -----------------------------

def read_hotspot_csv(path: Path):
    rows = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get("step") is None:
                # allow files that don't record step (legacy)
                step = 0
            else:
                try:
                    step = int(row["step"])
                except Exception:
                    step = 0
            rows.append({
                "step": step,
                "resname": row["resname"].strip().upper(),
                "resid":   int(row["resid"]),
                "x": float(row["x_nm"]),
                "y": float(row["y_nm"]),
                "z": float(row["z_nm"]),
            })
    return rows

def group_time_averages(rows):
    # average each (resname,resid) trajectory → one time-avg point
    acc = {}
    cnt = {}
    for r in rows:
        key = (r["resname"], r["resid"])
        if key not in acc:
            acc[key] = np.zeros(3, dtype=float)
            cnt[key] = 0
        acc[key] += np.array([r["x"], r["y"], r["z"]], dtype=float)
        cnt[key] += 1
    out = defaultdict(list)  # resname -> [point_nm]
    for (resname, resid), v in acc.items():
        if cnt[(resname, resid)] > 0:
            out[resname].append(v / cnt[(resname, resid)])
    return out  # dict: name -> [points_nm]

# -----------------------------
# Greedy clustering (dependency-free core)
# -----------------------------

def greedy_cluster(points_nm, radius_nm):
    """Return clusters as list of (centroid_nm, member_indices)."""
    if not points_nm:
        return []
    pts = np.array(points_nm, dtype=float)  # (N,3)
    used = np.zeros(len(pts), dtype=bool)
    clusters = []
    r2 = float(radius_nm) ** 2

    for i in range(len(pts)):
        if used[i]:
            continue
        # seed a new cluster
        center = pts[i].copy()
        members = [i]
        used[i] = True

        # single-pass greedy merge
        changed = True
        while changed:
            changed = False
            for j in range(len(pts)):
                if used[j]:
                    continue
                # distance to current center
                d2 = np.sum((pts[j] - center) ** 2)
                if d2 <= r2:
                    members.append(j)
                    used[j] = True
                    # update centroid
                    center = (center * (len(members) - 1) + pts[j]) / len(members)
                    changed = True

        clusters.append((center, members))
    # sort by size (desc)
    clusters.sort(key=lambda c: len(c[1]), reverse=True)
    return clusters

# -----------------------------
# PDB writing
# -----------------------------

PROBE_TO_CHAIN = ["I","A","M","C","P","H","U","V","W","X","Y","Z"]  # enough variety

def _fmt_pdb_atom(serial, name, resn, chain, resid, xA, yA, zA, bfac=1.0, occ=1.0):
    # name: 4 chars right/aligned (we use 'CEN ')
    return (f"HETATM{serial:5d} {name:<4s}{resn:>3s} {chain:1s}"
            f"{resid:>4d}    {xA:8.3f}{yA:8.3f}{zA:8.3f}{occ:6.2f}{bfac:6.2f}          C  \n")

def write_hotspot_pdb_trajframe(out_path: Path, per_probe_clusters, title_remark):
    """Write centroids as standalone PDB in the *trajectory* frame (Å)."""
    serial = 1
    with open(out_path, "w") as f:
        f.write("REMARK 500 HOTSPOT CENTROIDS (trajectory frame)\n")
        f.write("REMARK 500 B-FACTOR FIELD = number of probe instances in cluster\n")
        if title_remark:
            f.write(f"REMARK 500 {title_remark}\n")
        chain_map = {}
        for idx, resn in enumerate(per_probe_clusters.keys()):
            chain_map[resn] = PROBE_TO_CHAIN[idx % len(PROBE_TO_CHAIN)]
        for resn, clusters in per_probe_clusters.items():
            chain = chain_map[resn]
            for rank, (cen_nm, members) in enumerate(clusters, start=1):
                xA, yA, zA = (cen_nm*10.0).tolist()  # nm → Å
                f.write(_fmt_pdb_atom(serial, "CEN", resn, chain, rank, xA, yA, zA,
                                      bfac=float(len(members)), occ=1.00))
                serial += 1
            f.write("TER\n")
        f.write("END\n")

# -----------------------------
# Alignment & merged PDB
# -----------------------------

def _parse_ca_coords(pdb_path: Path):
    """Return Nx3 array of CA atom coordinates (Å) in file order."""
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    coords.append([x, y, z])
                except Exception:
                    pass
    return np.array(coords, dtype=float)

def _kabsch(P, Q):
    """Find rotation R, translation t such that P*R + t ≈ Q (both Nx3)."""
    # centroids
    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    P0 = P - cP
    Q0 = Q - cQ
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # ensure right-handed
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = cQ - cP @ R
    return R, t

def write_hotspot_pdb_on_reference(prefix: Path,
                                   per_probe_clusters,
                                   trajframe_pdb: Path,
                                   reference_pdb: Path):
    """
    Align trajframe Cα → reference Cα and merge centroids onto reference.
    Outputs <prefix>_on_<reference-stem>.pdb
    """
    # CA sets (Å)
    P = _parse_ca_coords(trajframe_pdb)
    Q = _parse_ca_coords(reference_pdb)
    n = min(len(P), len(Q))
    if n < 5:
        print("⚠️  Not enough matching CA atoms to align; skip merged PDB.")
        return None
    R, t = _kabsch(P[:n], Q[:n])

    # transform centroids (Å)
    per_probe_A = {}
    for resn, clusters in per_probe_clusters.items():
        out = []
        for cen_nm, members in clusters:
            vA = cen_nm * 10.0  # nm→Å
            vA = vA @ R + t
            out.append((vA, members))
        per_probe_A[resn] = out

    # read reference PDB to copy atoms
    ref_lines = []
    with open(reference_pdb, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TER", "MODEL", "ENDMDL", "REMARK", "HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA")):
                ref_lines.append(line)
    # determine next serial
    max_serial = 0
    for line in ref_lines:
        if line.startswith(("ATOM","HETATM")):
            try:
                s = int(line[6:11])
                if s > max_serial: max_serial = s
            except Exception:
                pass
    serial = max_serial + 1

    out_path = prefix.with_name(f"{prefix.name}_on_{Path(reference_pdb).stem}.pdb")
    with open(out_path, "w") as f:
        f.write(f"REMARK 500 HOTSPOT CENTROIDS merged onto reference: {Path(reference_pdb).name}\n")
        f.write("REMARK 500 B-FACTOR FIELD = number of probe instances in cluster\n")
        # write reference first
        for line in ref_lines:
            if not line.startswith("END"):  # avoid duplicate END
                f.write(line)
        # append centroids
        chain_map = {}
        for idx, resn in enumerate(per_probe_A.keys()):
            chain_map[resn] = PROBE_TO_CHAIN[idx % len(PROBE_TO_CHAIN)]
        for resn, items in per_probe_A.items():
            chain = chain_map[resn]
            for rank, (vA, members) in enumerate(items, start=1):
                xA, yA, zA = vA.tolist()
                f.write(_fmt_pdb_atom(serial, "CEN", resn, chain, rank, xA, yA, zA,
                                      bfac=float(len(members)), occ=1.00))
                serial += 1
            f.write("TER\n")
        f.write("END\n")
    return out_path

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Greedy cluster MixMD probe centroids and write hotspots (CSV/PDB)."
    )
    ap.add_argument("--input", required=True, help="mixmd_hotspots*.csv recorded by 5_md_simulation.py")
    ap.add_argument("--radius-nm", type=float, default=0.40, help="merge radius (nm) for greedy clustering")
    ap.add_argument("--min-members", type=int, default=3, help="minimum members to keep a cluster")
    ap.add_argument("--top-k", type=int, default=3, help="top clusters per probe type to write")

    # Alignment / merging options
    ap.add_argument("--source-pdb", type=str, default=None,
                    help="PDB in the trajectory frame (auto: final_structure*_no_ligand.pdb). Used for alignment.")
    ap.add_argument("--reference-pdb", type=str, default=None,
                    help="Reference receptor PDB to merge onto (auto: combined_receptor_ligand*_no_ligand.pdb, "
                         "otherwise receptor_cleaned.pdb).")
    ap.add_argument("--no-merge", action="store_true",
                    help="Do not write the merged-on-reference PDB; only write the trajectory-frame PDB.")

    args = ap.parse_args()
    in_csv = Path(args.input)

    rows = read_hotspot_csv(in_csv)
    if not rows:
        raise SystemExit(f"No rows parsed from {in_csv}")

    # time-averaged points per instance (per probe type)
    per_probe_points = group_time_averages(rows)

    # cluster per probe & collect top-K (by size)
    per_probe_clusters = {}
    vina_boxes = []  # for optional text summary
    for resn, pts in per_probe_points.items():
        clusters = greedy_cluster(pts, args.radius_nm)
        # filter by min-members
        clusters = [(c, m) for (c, m) in clusters if len(m) >= args.min_members]
        if not clusters:
            continue
        # keep top-K
        clusters = clusters[: args.top_k]
        per_probe_clusters[resn] = clusters
        # quick Vina boxes (Å): center, 22 Å default cube
        for rank, (cen_nm, members) in enumerate(clusters, start=1):
            xA, yA, zA = (cen_nm * 10.0).tolist()
            vina_boxes.append((resn, rank, xA, yA, zA, len(members)))

    if not per_probe_clusters:
        raise SystemExit("No clusters survived filtering; try a larger --radius-nm or lower --min-members.")

    # filenames
    prefix = Path(f"hotspots_{in_csv.stem}_r{args.radius_nm:.2f}nm")

    # 1) Trajectory-frame PDB (standalone)
    write_hotspot_pdb_trajframe(
        prefix.with_name(f"{prefix.name}_trajframe.pdb"),
        per_probe_clusters,
        title_remark=f"radius={args.radius_nm:.2f} nm; min_members={args.min_members}; top_k={args.top_k}"
    )

    # 2) (optional) merged onto reference receptor PDB
    if not args.no_merge:
        # infer defaults
        src = Path(args.source_pdb) if args.source_pdb else None
        ref = Path(args.reference_pdb) if args.reference_pdb else None

        if src is None:
            for cand in ["final_structure_no_ligand.pdb", "final_structure.pdb",
                         "npt_equilibrated_no_ligand.pdb", "nvt_equilibrated_no_ligand.pdb"]:
                if Path(cand).exists():
                    src = Path(cand); break
        if ref is None:
            for cand in ["combined_receptor_ligand_no_ligand.pdb",
                         "combined_receptor_ligand.pdb",
                         "receptor_cleaned.pdb"]:
                if Path(cand).exists():
                    ref = Path(cand); break

        if src and ref and src.exists() and ref.exists():
            merged = write_hotspot_pdb_on_reference(prefix, per_probe_clusters, src, ref)
            if merged:
                print(f"✅ Wrote merged hotspots PDB → {merged}")
            else:
                print("⚠️  Merged PDB not written (alignment failed); see warnings above.")
        else:
            print(f"⚠️  Missing source/ref PDB(s) for alignment (src={src}, ref={ref}); only trajframe PDB written.")

    # 3) minimal AutoDock‑Vina boxes summary
    with open(prefix.with_suffix(".vina.txt"), "w") as f:
        f.write("# center_x  center_y  center_z  size_x  size_y  size_z  ;  resname  rank  members\n")
        for resn, rank, xA, yA, zA, n in vina_boxes:
            f.write(f"{xA:8.3f}  {yA:8.3f}  {zA:8.3f}   22.0  22.0  22.0  ;  {resn}  {rank}  {n}\n")

    # 4) CSV summary (same as before)
    out_csv = prefix.with_suffix(".csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["probe", "rank", "members", "x_A", "y_A", "z_A"])
        for resn, clusters in per_probe_clusters.items():
            for rank, (cen_nm, members) in enumerate(clusters, start=1):
                xA, yA, zA = (cen_nm * 10.0).tolist()
                w.writerow([resn, rank, len(members), f"{xA:.3f}", f"{yA:.3f}", f"{zA:.3f}"])

    print(f"✅ Hotspots CSV   → {out_csv}")
    print(f"✅ Hotspots PDB   → {prefix.name}_trajframe.pdb")
    if not args.no_merge and Path(f"{prefix.name}_on_{Path(ref).stem}.pdb").exists():
        print(f"✅ Hotspots PDB   → {prefix.name}_on_{Path(ref).stem}.pdb (aligned & merged)")

if __name__ == "__main__":
    main()
