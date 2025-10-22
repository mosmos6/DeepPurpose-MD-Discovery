# 5c_mixmd_hotspot_cluster.py
# Greedy clustering of MixMD probe centroids with optional PBC wrapping and Kabsch alignment
# Updated: 2025-10-22

import argparse, csv
from pathlib import Path
from collections import defaultdict

import numpy as np

# -----------------------------
# I/O helpers
# -----------------------------

def read_hotspot_csv(path: Path):
    """
    Read hotspots CSV written by 5 (ProbeHotspotReporter).
    - Skips comment lines (starting with '#').
    - Accepts header aliases: x_nm|x, y_nm|y, z_nm|z, resid|res_id|id, resname|probe|name
    - 'step' is optional.
    Returns: list of dicts with keys: step,resname,resid,x,y,z (nm floats).
    """
    rows = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        clean = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    if not clean:
        return rows

    rdr = csv.DictReader(clean)
    for rec in rdr:
        lower = { (k or "").strip().lower(): (v or "").strip()
                  for k, v in rec.items() }

        def pick(*names, default=None):
            for n in names:
                if n in lower and lower[n] != "":
                    return lower[n]
            return default

        resname = (pick("resname", "probe", "name", default="")).upper()
        x_s = pick("x_nm", "x"); y_s = pick("y_nm", "y"); z_s = pick("z_nm", "z")
        resid_s = pick("resid", "res_id", "id", default="0")
        step_s  = pick("step", default="0")
        if not (resname and x_s is not None and y_s is not None and z_s is not None):
            continue
        try:
            x = float(x_s); y = float(y_s); z = float(z_s)
            resid = int(float(resid_s)); step = int(float(step_s))
        except Exception:
            continue
        rows.append({"step": step, "resname": resname, "resid": resid, "x": x, "y": y, "z": z})
    return rows


def group_time_averages(rows):
    """Average each (resname,resid) trajectory → one point (nm)."""
    acc, cnt = {}, {}
    for r in rows:
        key = (r["resname"], r["resid"])
        if key not in acc:
            acc[key] = np.zeros(3, float); cnt[key] = 0
        acc[key] += np.array([r["x"], r["y"], r["z"]], float)
        cnt[key] += 1
    out = defaultdict(list)  # resname -> [point_nm]
    for (resname, resid), v in acc.items():
        if cnt[(resname, resid)]:
            out[resname].append(v / cnt[(resname, resid)])
    return out

# -----------------------------
# PBC wrapping near protein center
# -----------------------------

def _protein_ca_center_nm(pdb_path: Path):
    """Return mean Cα position (nm) from a PDB."""
    ca = []
    with open(pdb_path, "r") as f:
        for ln in f:
            if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                    ca.append((x, y, z))
                except Exception:
                    pass
    if not ca:
        return np.zeros(3, float)
    arr = np.asarray(ca, float) / 10.0  # Å → nm
    return arr.mean(axis=0)

def _wrap_to_near(v_nm, ref_nm, L_nm):
    """Wrap coordinate v_nm (3,) to the image nearest ref_nm given edge L_nm."""
    return ref_nm + ((v_nm - ref_nm + 0.5*L_nm) % L_nm) - 0.5*L_nm

def wrap_rows_near_protein(rows, L_nm: float, src_pdb_for_center: Path | None):
    """In-place wrap of each row's (x,y,z) so it lies near the protein center."""
    if src_pdb_for_center and src_pdb_for_center.exists():
        ctr = _protein_ca_center_nm(src_pdb_for_center)
    else:
        # last resort: use mean of recorded points
        pts = np.array([[r["x"], r["y"], r["z"]] for r in rows], float)
        ctr = pts.mean(axis=0) if len(pts) else np.zeros(3, float)

    for r in rows:
        v = np.array([r["x"], r["y"], r["z"]], float)
        w = _wrap_to_near(v, ctr, float(L_nm))
        r["x"], r["y"], r["z"] = map(float, w.tolist())

# -----------------------------
# Greedy clustering
# -----------------------------

def greedy_cluster(points_nm, radius_nm):
    """Return clusters as list of (centroid_nm, member_indices)."""
    if not points_nm:
        return []
    pts = np.array(points_nm, float)
    used = np.zeros(len(pts), bool)
    clusters, r2 = [], float(radius_nm) ** 2

    for i in range(len(pts)):
        if used[i]:
            continue
        center = pts[i].copy()
        members = [i]
        used[i] = True
        changed = True
        while changed:
            changed = False
            for j in range(len(pts)):
                if used[j]:
                    continue
                d2 = float(np.sum((pts[j] - center) ** 2))
                if d2 <= r2:
                    members.append(j); used[j] = True
                    center = (center * (len(members) - 1) + pts[j]) / len(members)
                    changed = True
        clusters.append((center, members))
    clusters.sort(key=lambda c: len(c[1]), reverse=True)
    return clusters

# -----------------------------
# PDB writing (correct columns)
# -----------------------------

PROBE_TO_CHAIN = ["I","A","M","C","P","H","U","V","W","X","Y","Z"]  # variety

def _fmt_pdb_atom(serial, name, resn, chain, resid, xA, yA, zA, bfac=1.0, occ=1.0, altloc=" "):
    # Columns: 1-6 rec, 7-11 serial, 13-16 name, 17 altLoc, 18-20 resName, 22 chain, 23-26 resSeq
    return (f"HETATM{serial:5d} {name:<4s}{altloc:1s}{resn:>3s} {chain:1s}"
            f"{resid:4d}    {xA:8.3f}{yA:8.3f}{zA:8.3f}{occ:6.2f}{bfac:6.2f}          C  \n")

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
                xA, yA, zA = (cen_nm * 10.0).tolist()
                f.write(_fmt_pdb_atom(serial, "CEN", resn, chain, rank, xA, yA, zA,
                                      bfac=float(len(members)), occ=1.00, altloc=" "))
                serial += 1
            f.write("TER\n")
        f.write("END\n")

# -----------------------------
# Alignment & merged PDB
# -----------------------------

def _parse_ca_indexed(pdb_path: Path):
    """Return {(chain, resSeq, iCode) -> [x,y,z]Å} and an ordered array for NN fallback."""
    ca_dict, ca_list = {}, []
    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):  continue
            if line[12:16].strip() != "CA":   continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):       continue
            try:
                chain = line[21].strip() or "_"
                resseq = int(line[22:26])
                icode  = line[26].strip() or ""
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except Exception:
                continue
            key = (chain, resseq, icode)
            v = np.array([x, y, z], float)
            ca_dict[key] = v; ca_list.append(v)
    arr = np.vstack(ca_list) if ca_list else np.zeros((0,3), float)
    return ca_dict, arr

def _pair_ca_by_ids(src_dict, ref_dict):
    keys = sorted(set(src_dict.keys()) & set(ref_dict.keys()),
                  key=lambda k: (k[0], k[1], k[2]))
    if not keys:
        return np.zeros((0,3), float), np.zeros((0,3), float)
    P = np.array([src_dict[k] for k in keys], float)
    Q = np.array([ref_dict[k] for k in keys], float)
    return P, Q

def _pair_ca_by_nearest(src_all, ref_all, max_dist_A=5.0):
    if src_all.size == 0 or ref_all.size == 0:
        return np.zeros((0,3), float), np.zeros((0,3), float)
    s2 = np.sum(src_all**2, axis=1)[:, None]
    r2 = np.sum(ref_all**2, axis=1)[None, :]
    D = s2 + r2 - 2.0 * (src_all @ ref_all.T)
    used_ref, pairs = set(), []
    for i in range(D.shape[0]):
        j = int(np.argmin(D[i])); d2 = float(D[i, j])
        if j not in used_ref and d2 <= max_dist_A*max_dist_A:
            pairs.append((i, j)); used_ref.add(j)
    if not pairs:
        return np.zeros((0,3), float), np.zeros((0,3), float)
    P = np.array([src_all[i] for i, _ in pairs], float)
    Q = np.array([ref_all[j] for _, j in pairs], float)
    return P, Q

def _kabsch(P, Q):
    cP, cQ = P.mean(axis=0), Q.mean(axis=0)
    P0, Q0 = P - cP, Q - cQ
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1; R = Vt.T @ U.T
    t = cQ - cP @ R
    return R, t

def _align_pairs_iterative(P, Q, prune_A=3.0, iters=3):
    if len(P) < 3 or len(Q) < 3:
        return np.eye(3), np.zeros(3), np.inf, 0
    mask = np.ones(len(P), bool)
    R = np.eye(3); t = np.zeros(3)
    for _ in range(iters):
        if mask.sum() < 3: break
        R, t = _kabsch(P[mask], Q[mask])
        P_aln = (P @ R) + t
        d = np.linalg.norm(P_aln - Q, axis=1)
        mask = d < prune_A
    rmsd = float(np.sqrt(np.mean(((P @ R + t) - Q)**2)))
    return R, t, rmsd, int(mask.sum())

def write_hotspot_pdb_on_reference(prefix: Path,
                                   per_probe_clusters,
                                   trajframe_pdb: Path,
                                   reference_pdb: Path):
    """Align trajframe Cα → reference Cα and merge centroids onto the reference (Å)."""
    src_dict, src_all = _parse_ca_indexed(trajframe_pdb)
    ref_dict, ref_all = _parse_ca_indexed(reference_pdb)
    P, Q = _pair_ca_by_ids(src_dict, ref_dict); mode = "id"
    if len(P) < 5:
        P, Q = _pair_ca_by_nearest(src_all, ref_all, max_dist_A=5.0); mode = "nn"
    if len(P) < 5:
        print(f"⚠️  Alignment failed: only {len(P)} CA pairs (mode={mode}). Skipping merged PDB.")
        return None
    pre = float(np.sqrt(np.mean((P - Q)**2)))
    R, t, post, used = _align_pairs_iterative(P, Q, prune_A=3.0, iters=3)
    print(f"[align] mode={mode} pairs={len(P)} used={used}  RMSD_before={pre:.2f} Å  RMSD_after={post:.2f} Å")

    # transform centroids (nm → Å) and apply R,t
    per_probe_A = {}
    for resn, clusters in per_probe_clusters.items():
        out = []
        for cen_nm, members in clusters:
            vA = (cen_nm * 10.0) @ R + t
            out.append((vA, members))
        per_probe_A[resn] = out

    # copy reference lines, find next serial
    ref_lines, max_serial = [], 0
    with open(reference_pdb, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TER", "MODEL", "ENDMDL",
                                "REMARK", "HEADER", "TITLE", "COMPND", "SOURCE",
                                "KEYWDS", "EXPDTA")):
                ref_lines.append(line)
                if line.startswith(("ATOM", "HETATM")):
                    try: max_serial = max(max_serial, int(line[6:11]))
                    except Exception: pass
    serial = max_serial + 1

    out_path = prefix.with_name(f"{prefix.name}_on_{Path(reference_pdb).stem}.pdb")
    chain_map = {}
    with open(out_path, "w") as f:
        f.write(f"REMARK 500 HOTSPOT CENTROIDS merged onto reference: {Path(reference_pdb).name}\n")
        f.write("REMARK 500 B-FACTOR FIELD = number of probe instances in cluster\n")
        for line in ref_lines:
            if not line.startswith("END"):
                f.write(line)
        for idx, resn in enumerate(per_probe_A.keys()):
            chain_map[resn] = PROBE_TO_CHAIN[idx % len(PROBE_TO_CHAIN)]
        for resn, items in per_probe_A.items():
            chain = chain_map[resn]
            for rank, (vA, members) in enumerate(items, start=1):
                xA, yA, zA = vA.tolist()
                f.write(_fmt_pdb_atom(serial, "CEN", resn, chain, rank, xA, yA, zA,
                                      bfac=float(len(members)), occ=1.00, altloc=" "))
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
    ap.add_argument("--radius-nm", type=float, default=0.40, help="merge radius (nm)")
    ap.add_argument("--min-members", type=int, default=3, help="minimum members to keep a cluster")
    ap.add_argument("--top-k", type=int, default=3, help="top clusters per probe to write")

    # Optional PBC wrapping (highly recommended for MixMD)
    ap.add_argument("--box-size-nm", type=float, default=None,
                    help="If set, wrap each recorded point to the image nearest the protein center.")
    ap.add_argument("--source-pdb", type=str, default=None,
                    help="Trajectory-frame PDB for protein center & alignment (auto: final_structure*_no_ligand.pdb).")
    ap.add_argument("--reference-pdb", type=str, default=None,
                    help="Reference receptor PDB for merged output (auto: combined_receptor*_no_ligand.pdb, else receptor_cleaned.pdb).")
    ap.add_argument("--no-merge", action="store_true",
                    help="Do not write the merged-on-reference PDB; only write the trajectory-frame PDB.")

    args = ap.parse_args()
    in_csv = Path(args.input)

    rows = read_hotspot_csv(in_csv)
    if not rows:
        raise SystemExit(f"No rows parsed from {in_csv}")

    # Optional: wrap to nearest protein image (fixes 'cluster in the sky')
    src_pdb = Path(args.source_pdb) if args.source_pdb else None
    if src_pdb is None:
        for cand in ["final_structure_no_ligand.pdb", "final_structure.pdb",
                     "npt_equilibrated_no_ligand.pdb", "nvt_equilibrated_no_ligand.pdb"]:
            if Path(cand).exists():
                src_pdb = Path(cand); break

    if args.box_size_nm is not None:
        wrap_rows_near_protein(rows, float(args.box_size_nm), src_pdb)

    # Time-averaged points per instance
    per_probe_points = group_time_averages(rows)

    # Cluster per probe
    per_probe_clusters = {}
    vina_boxes = []
    for resn, pts in per_probe_points.items():
        clusters = greedy_cluster(pts, args.radius_nm)
        clusters = [(c, m) for (c, m) in clusters if len(m) >= args.min_members]
        if not clusters:
            continue
        clusters = clusters[: args.top_k]
        per_probe_clusters[resn] = clusters
        for rank, (cen_nm, members) in enumerate(clusters, start=1):
            xA, yA, zA = (cen_nm * 10.0).tolist()
            vina_boxes.append((resn, rank, xA, yA, zA, len(members)))

    if not per_probe_clusters:
        raise SystemExit("No clusters survived filtering; try a larger --radius-nm or lower --min-members.")

    # Filenames
    prefix = Path(f"hotspots_{in_csv.stem}_r{args.radius_nm:.2f}nm")

    # 1) Trajectory-frame PDB
    write_hotspot_pdb_trajframe(
        prefix.with_name(f"{prefix.name}_trajframe.pdb"),
        per_probe_clusters,
        title_remark=f"radius={args.radius_nm:.2f} nm; min_members={args.min_members}; top_k={args.top_k}"
    )

    # 2) Merged-on-reference PDB (aligned)
    if not args.no_merge:
        ref_pdb = Path(args.reference_pdb) if args.reference_pdb else None
        if ref_pdb is None:
            for cand in ["combined_receptor_ligand_no_ligand.pdb",
                         "combined_receptor_ligand.pdb",
                         "receptor_cleaned.pdb"]:
                if Path(cand).exists():
                    ref_pdb = Path(cand); break
        if src_pdb and ref_pdb and src_pdb.exists() and ref_pdb.exists():
            merged = write_hotspot_pdb_on_reference(prefix, per_probe_clusters, src_pdb, ref_pdb)
            if merged:
                print(f"✅ Wrote merged hotspots PDB → {merged}")
            else:
                print("⚠️  Merged PDB not written (alignment failed); see warnings above.")
        else:
            print(f"⚠️  Missing source/ref PDB(s) for alignment (src={src_pdb}, ref={ref_pdb}); only trajframe PDB written.")

    # 3) Vina boxes (Å)
    with open(prefix.with_suffix(".vina.txt"), "w") as f:
        f.write("# center_x  center_y  center_z  size_x  size_y  size_z  ;  resname  rank  members\n")
        for resn, rank, xA, yA, zA, n in vina_boxes:
            f.write(f"{xA:8.3f}  {yA:8.3f}  {zA:8.3f}   22.0  22.0  22.0  ;  {resn}  {rank}  {n}\n")

    # 4) CSV summary
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

if __name__ == "__main__":
    main()
