# 2_prepare_receptor.py
# Author: Iori Mochizuki
# Updated: 2025-09-11
# Description: Prepare protein (by PDB ID) or RNA (SimRNA PDB) for docking and MD

import os, sys, subprocess, argparse

parser = argparse.ArgumentParser(description="Prepare receptor (protein or RNA) for docking and MD.")
parser.add_argument("--pdb-id", type=str, default=None, help="PDB ID to download from RCSB (protein mode)")
parser.add_argument("--rna", action="store_true", help="Enable RNA mode: use uploaded PDB file (SimRNA output)")
parser.add_argument("--input-pdb", type=str, default="receptor.pdb", help="Input RNA PDB file if --rna (default: receptor.pdb)")
parser.add_argument("--skip-fix", action="store_true", help="Skip PDBFixer step (protein mode only)")
parser.add_argument("--strict-protein", action="store_true", help="Strip non-protein residues using MDTraj (protein mode only)")

# NEW: resolve altloc-style duplicated atom sites (multi-chain safe). Off by default.
parser.add_argument("--resolve-altdups", action="store_true",
                    help="Resolve duplicated atom sites across all chains by keeping the highest-occupancy record per site (handles blank altLoc).")

args = parser.parse_args()

receptor_clean = "receptor_clean.pdb"
receptor_for_centroid = "receptor_for_centroid.pdb"


def resolve_multichain_alt_dups(input_pdb: str, output_pdb: str, keep_hetatm: bool = True, keep_water: bool = False):
    """
    Collapse duplicated atom sites across ALL chains in a PDB.
    - Keeps highest-occupancy record per (chain, resseq, icode, atomname).
    - If occupancy ties, prefer blank altLoc; otherwise keep first encountered.
    - Optionally drop HETATM and/or HOH.
    """
    def is_atom(line):
        return line.startswith("ATOM") or line.startswith("HETATM")

    def is_water(line):
        return line[17:20].strip() == "HOH"

    def is_het(line):
        return line.startswith("HETATM")

    best = {}       # key -> (occ, prefer_blank_altloc, line)
    other = []      # non-ATOM/HETATM lines (HEADER/REMARK/TER/END, etc.)

    with open(input_pdb, "r") as fh:
        for line in fh:
            if not is_atom(line):
                other.append(line)
                continue

            if is_water(line) and not keep_water:
                continue
            if is_het(line) and not keep_hetatm:
                continue

            chain  = line[21]
            resseq = line[22:26]
            icode  = line[26]
            aname  = line[12:16]
            altloc = line[16]  # may be blank
            try:
                occ = float(line[54:60])
            except ValueError:
                occ = 0.0

            key = (chain, resseq, icode, aname)
            cand = (occ, (altloc == " "), line)
            if key not in best:
                best[key] = cand
            else:
                occ0, prefer_blank0, _line0 = best[key]
                # higher occupancy wins; if tie, prefer blank altLoc
                if cand[0] > occ0 or (cand[0] == occ0 and cand[1] and not prefer_blank0):
                    best[key] = cand

    chosen = []
    for _, (_, _, line) in best.items():
        # normalize: blank out altLoc column (col 17)
        fixed = line[:16] + " " + line[17:]
        chosen.append(fixed)

    safe_header = [
        ln for ln in other
        if not (ln.startswith("CONECT") or ln.startswith("END") or ln.startswith("ENDMDL") or ln.startswith("MASTER"))
    ]

    with open(output_pdb, "w") as out:
        # 1) safe headers / remarks only
        for line in safe_header:
            out.write(line)
        # 2) atom records
        for line in chosen:
            out.write(line)
        # 3) write a single clean END at the very end
        out.write("END\n")

    print(f"🔧 Alt/dup resolver → {output_pdb} (kept {len(chosen)} atom records)")


if args.rna:
    # ---- RNA MODE ----
    print("🧬 [RNA MODE] Using uploaded PDB:", args.input_pdb)
    # Copy the input RNA PDB to standard names (no changes to behavior)
    subprocess.run(["cp", args.input_pdb, receptor_clean])
    subprocess.run(["cp", args.input_pdb, receptor_for_centroid])
    receptor_for_meeko = receptor_clean

else:
    # ---- PROTEIN MODE ----
    if not args.pdb_id:
        print("❌ Error: --pdb-id must be specified unless --rna is set.")
        sys.exit(1)

    receptor_raw = "receptor.pdb"
    receptor_fixed = "receptor_fixed.pdb"

    print(f"📦 Downloading {args.pdb_id} from RCSB...")
    #subprocess.run(["wget", f"https://files.rcsb.org/download/{args.pdb_id}.pdb", "-O", receptor_raw], check=False)

    # 1) Make two branches from the raw file, removing HOH in both
    #    - receptor_for_centroid.pdb : keep HETATM (for centroid logic)
    #    - receptor_clean_base.pdb   : will become docking/MD-ready (no HET)
    print("🧹 Removing water molecules...")
    subprocess.run(["obabel", receptor_raw, "-O", receptor_for_centroid, "--delete", "HOH"], check=False)
    subprocess.run(["obabel", receptor_raw, "-O", "receptor_clean_base.pdb", "--delete", "HOH"], check=False)

    # 2) Resolve multi-chain duplicated atom sites / altloc-style conformers (optional flag)
    if args.resolve_altdups:
        print("🧼 Resolving duplicates (centroid branch, KEEP HETATM)…")
        resolve_multichain_alt_dups(
            input_pdb=receptor_for_centroid,
            output_pdb=receptor_for_centroid,   # overwrite in place
            keep_hetatm=True,
            keep_water=False
        )

        print("🧼 Resolving duplicates (clean branch, DROP HETATM)…")
        resolve_multichain_alt_dups(
            input_pdb="receptor_clean_base.pdb",
            output_pdb="receptor_clean_nohet_noalt.pdb",
            keep_hetatm=False,
            keep_water=False
        )
        # adopt deduped file as the clean branch
        receptor_clean_path = "receptor_clean_nohet_noalt.pdb"
    else:
        # no resolver: pass through the HOH-stripped base as "clean" (still has HET)
        receptor_clean_path = "receptor_clean_base.pdb"

    # 3) Optional: strict-protein slice (applies to the clean branch only)
    if args.strict_protein:
        print("🔬 Removing non-protein residues using MDTraj (clean branch)…")
        import mdtraj as md
        traj = md.load_pdb(receptor_clean_path)
        filtered = traj.atom_slice(traj.top.select("protein"))
        filtered.save_pdb("receptor_clean.pdb")
        print("✅ Saved as receptor_clean.pdb")
    else:
        # rename the current clean to expected name for downstream compatibility
        if receptor_clean_path != "receptor_clean.pdb":
            os.replace(receptor_clean_path, "receptor_clean.pdb")

    # Ensure we point to the canonical output name
    receptor_clean = "receptor_clean.pdb"

    # 4) Optional PDBFixer (on the clean branch). Centroid branch is untouched (by design).
    if args.skip_fix:
        print("⚠️ Skipping PDBFixer cleanup (as requested via --skip-fix).")
        receptor_for_meeko = receptor_clean
    else:
        print("🧬 Running PDBFixer on receptor_clean.pdb…")
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=receptor_clean)
        fixer.findMissingResidues()
        # Do NOT rebuild loops—keep native coordinates for high-res crystals
        fixer.missingResidues = {}

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)

        with open(receptor_fixed, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        print(f"✅ PDBFixer output → {receptor_fixed}")
        receptor_for_meeko = receptor_fixed

# --- Meeko: Always run for both protein and RNA ---
print("🛠️ Running Meeko on:", receptor_for_meeko)
meeko_command = (
    "python3 /usr/local/envs/deeppurpose-md-env/lib/python3.11/site-packages/meeko/cli/mk_prepare_receptor.py "
    f"-i {receptor_for_meeko} -o receptor -p -j -v "
    "--box_size 20 20 20 --box_center 0 0 0 --allow_bad_res"
)
os.system(meeko_command)
print("✅ Receptor PDBQT prepared.")
