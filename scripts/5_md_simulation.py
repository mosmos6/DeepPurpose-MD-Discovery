#!/usr/bin/env python3
# 5_md_simulation.py
# Author: Iori Mochizuki (+ collaborator)
# Updated: 2025-10-23
#
# OpenMM 2 ns MD (protein/RNA/ligand) with optional MixMD-from-Packmol import.
# - Single- or multi-ligand MD: auto-detects ligand*.sdf (or use --input-ligand).
# - MixMD-from-Packmol: places probes from build/<receptor>_mixmd_placements.csv,
#   renames residues (IPA/ACN/IMD/ACM/PHO/HAC), and logs probe hotspots each stride.
# - Hotspot CSV logs **beyond step 0** (robust scheduling + COM-centering).
# - Writes PDB **and** PDBx/mmCIF (for big boxes where PDB serials overflow).
#
# Example:
#   conda run -n deeppurpose-md-env python scripts/5_md_simulation.py \
#       --mixmd-from-packmol --mixmd-box-size-nm 8.0 --hotspot-stride 200 \
#       --seed 111 --n-steps 100000
#
# Multi-ligand:
#   Place SDFs as: ligand.sdf, ligand_1.sdf, ligand_2.sdf or in ligands/*.sdf
#   and run **without** --mixmd-from-packmol:
#   conda run -n deeppurpose-md-env python scripts/5_md_simulation.py --n-steps 100000

import argparse, csv, glob
from pathlib import Path
from sys import stdout

import openmm as mm
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---------------------------
# Probe library (MixMD)
# ---------------------------
PROBE_LIBRARY = {
    # canonical 3-letter → SMILES
    "IPA":  "CC(C)O",           # isopropanol
    "ACN":  "CC#N",             # acetonitrile
    "IMD":  "c1[nH]cnc1",       # imidazole (neutral, N1–H)
    "ACM":  "CC(=O)N",          # acetamide
    "PHO":  "c1ccc(cc1)O",      # phenol
    "HAC":  "CC(=O)O",          # acetic acid
    # aliases accepted from CLI for convenience
    "ACEA": "CC(=O)N",          # acetamide (alias of ACM)
    "PHOL": "c1ccc(cc1)O",      # phenol (alias of PHO)
    "ACOH": "CC(=O)O",          # acetic acid (alias of HAC)
}
PROBE_NAME_SET = set(PROBE_LIBRARY.keys())

def _default_packmol_csv(receptor_path: Path) -> Path:
    root = f"{receptor_path.stem}_mixmd"
    return Path("build") / f"{root}_placements.csv"

# ---------------------------
# Hotspot reporter (robust)
# ---------------------------
class ProbeHotspotReporter:
    """
    Write probe residue centroids (nm) every 'reportInterval' steps to CSV.
    Columns: step,resname,resid,x_nm,y_nm,z_nm

    We log in a receptor-centered frame to remove drift:
        x' = x - COM_protein(t) + COM_protein(0)
    COM uses protein Cα atoms if available, else protein heavy atoms.
    """
    def __init__(self, file_path, reportInterval, topology, probe_resnames=None, center_mode="protein-com"):
        self.file = open(file_path, "w")
        self.reportInterval = int(reportInterval)
        self.center_mode = center_mode
        self._have_ref_com = False
        self._ref_com = (0.0, 0.0, 0.0)

        # set of probe names present (auto) ∪ user list (if provided)
        present = {(res.name or "").upper().strip() for res in topology.residues()}
        auto = present & PROBE_NAME_SET
        if probe_resnames:
            user = {s.upper().strip() for s in probe_resnames}
            self.probe_res = auto | (user & PROBE_NAME_SET)
        else:
            self.probe_res = auto

        # groups to average per probe residue
        groups, resid_map = [], []
        protein_resname_set = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        self._protein_anchor_indices = []

        for res in topology.residues():
            rn = (res.name or "").upper().strip()
            # probe residues
            if rn in self.probe_res:
                idxs = [a.index for a in res.atoms()]
                if idxs:
                    groups.append(idxs)
                    resid_map.append((rn, int(res.id) if res.id is not None else 0))
            # protein Cα as anchors
            if rn in protein_resname_set:
                for a in res.atoms():
                    if (a.name or "").upper() == "CA":
                        self._protein_anchor_indices.append(a.index)

        # fallback: protein heavy atoms if no CA
        if not self._protein_anchor_indices:
            for res in topology.residues():
                rn = (res.name or "").upper().strip()
                if rn in protein_resname_set:
                    for a in res.atoms():
                        if a.element is None or a.element.symbol == "H":
                            continue
                        self._protein_anchor_indices.append(a.index)

        self.groups = groups
        self.resid_map = resid_map

        # header + diagnostic
        self.file.write("# tracking_resnames=" + ";".join(sorted(self.probe_res)) + "\n")
        self.file.write("# center_mode=" + str(self.center_mode) + "\n")
        self.file.write("step,resname,resid,x_nm,y_nm,z_nm\n")
        self.file.flush()

    def describeNextReport(self, simulation):
        # Return tuple: (interval, positions, velocities, forces, energies, enforcePeriodicBox)
        return (self.reportInterval, True, False, False, False, True)

    @staticmethod
    def _com_of_indices(positions, idxs):
        sx = sy = sz = 0.0
        n = 0
        for i in idxs:
            v = positions[i]
            sx += float(v.x); sy += float(v.y); sz += float(v.z)
            n += 1
        if n == 0:
            return (0.0, 0.0, 0.0)
        return (sx/n, sy/n, sz/n)

    def _current_step(self, simulation, state):
        # robust step retrieval across OpenMM versions
        try:
            return int(simulation.currentStep)          # OpenMM ≥ 7.5
        except Exception:
            try:
                t  = state.getTime()                    # Quantity
                dt = simulation.integrator.getStepSize()
                return int(round((t/dt)))
            except Exception:
                return -1

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=False)  # list of Vec3 (nm)
        step = self._current_step(simulation, state)

        # receptor-centered shift
        if self.center_mode == "protein-com" and self._protein_anchor_indices:
            com_now = self._com_of_indices(pos, self._protein_anchor_indices)
            if not self._have_ref_com:
                self._ref_com = com_now
                self._have_ref_com = True
            dx = self._ref_com[0] - com_now[0]
            dy = self._ref_com[1] - com_now[1]
            dz = self._ref_com[2] - com_now[2]
        else:
            dx = dy = dz = 0.0

        # write all probe centroids at this step
        for (idxs, (resname, resid)) in zip(self.groups, self.resid_map):
            sx = sy = sz = 0.0
            n = 0
            for i in idxs:
                v = pos[i]
                sx += float(v.x); sy += float(v.y); sz += float(v.z)
                n += 1
            if n > 0:
                self.file.write(f"{step},{resname},{resid},{sx/n+dx:.5f},{sy/n+dy:.5f},{sz/n+dz:.5f}\n")
        self.file.flush()

    def __del__(self):
        try:
            self.file.close()
        except Exception:
            pass

# ---------------------------
# MixMD placement helpers
# ---------------------------
def _read_centroids_csv(csv_path: Path):
    rows = []
    with open(csv_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({
                "resname": row["resname"].strip().upper(),
                "resid":   int(row["resid"]),
                "x_nm":    float(row["x_nm"]),
                "y_nm":    float(row["y_nm"]),
                "z_nm":    float(row["z_nm"]),
            })
    return rows

def _openff_conf_to_angstrom_list(off_mol: Molecule):
    """Return list[(x,y,z)_Å] from an OFF Molecule conformer using pure Python."""
    conf = off_mol.conformers[0]  # quantity with units
    try:
        arr = conf.to("angstrom").m  # Pint path
        return [(float(r[0]), float(r[1]), float(r[2])) for r in arr]
    except Exception:
        q = to_openmm(conf)  # OpenMM quantity
        out = []
        for v in q.value_in_unit(angstrom):
            out.append((float(v[0]), float(v[1]), float(v[2])))
        return out

def _positions_for_centroid_nm(off_mol: Molecule, centroid_nm_tuple):
    """Build OpenMM positions (list[Vec3]*nanometer) centered at given (x,y,z)_nm."""
    coords_A = _openff_conf_to_angstrom_list(off_mol)
    n = len(coords_A)
    if n == 0:
        raise ValueError("OFF molecule has no coordinates.")
    sx = sy = sz = 0.0
    for xA, yA, zA in coords_A:
        sx += xA; sy += yA; sz += zA
    cx = sx / n; cy = sy / n; cz = sz / n
    tx, ty, tz = (float(centroid_nm_tuple[0]),
                  float(centroid_nm_tuple[1]),
                  float(centroid_nm_tuple[2]))
    placed = []
    for xA, yA, zA in coords_A:
        xnm = (xA - cx) * 0.1 + tx  # Å→nm + translate
        ynm = (yA - cy) * 0.1 + ty
        znm = (zA - cz) * 0.1 + tz
        placed.append(mm.Vec3(xnm, ynm, znm))
    return [p * nanometer for p in placed]

def _add_probes_from_packmol(modeller,
                             forcefield,
                             receptor_path: Path,
                             sim_box_nm: float,
                             edge_margin_nm: float,
                             placements_csv: Path | None = None,
                             resname_list: list[str] | None = None) -> int:
    """
    Add probes at Packmol centroids cropped to the intended MD box.
    Renames newly-added residues to the probe resname for downstream reporting.
    """
    if placements_csv is None:
        placements_csv = _default_packmol_csv(receptor_path)
    print(f"ℹ️  MixMD inputs: CSV={placements_csv}")
    if not placements_csv.exists():
        raise FileNotFoundError(f"Placements CSV not found: {placements_csv}")

    half = 0.5 * float(sim_box_nm)
    keep_lo = - (half - float(edge_margin_nm))
    keep_hi = + (half - float(edge_margin_nm))

    whitelist = {r.upper().strip() for r in resname_list} if resname_list else None

    rows = []
    for rec in _read_centroids_csv(placements_csv):
        resname = rec["resname"]
        if resname not in PROBE_LIBRARY:
            continue
        if whitelist and resname not in whitelist:
            continue
        x, y, z = rec["x_nm"], rec["y_nm"], rec["z_nm"]
        if (keep_lo <= x <= keep_hi) and (keep_lo <= y <= keep_hi) and (keep_lo <= z <= keep_hi):
            rows.append((resname, rec["resid"], (x, y, z)))

    if not rows:
        print("⚠️  No probe centroids inside the target box; continuing with receptor-only.")
        return 0

    used_res = sorted({res for (res, _, _) in rows})
    off_mols = []
    for res in used_res:
        smi = PROBE_LIBRARY[res]
        off = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = res
        off_mols.append(off)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    top_cache = {m.name: OFFTopology.from_molecules([m]).to_openmm() for m in off_mols}
    mol_cache = {m.name: m for m in off_mols}

    def _res_list():
        return list(modeller.topology.residues())

    added = 0
    for (resname, pack_resid, cen_nm) in rows:
        top = top_cache[resname]
        off = mol_cache[resname]
        pos = _positions_for_centroid_nm(off, cen_nm)

        before_res = len(_res_list())
        modeller.add(top, pos)
        after = _res_list()

        new_res = after[before_res:]
        for r in new_res:
            try:
                r.name = resname
            except Exception:
                pass
            try:
                r.id = str(pack_resid)
            except Exception:
                pass
        added += 1

    print(f"✅ Placed {added} probe molecules (cropped to {sim_box_nm:.1f} nm box).")
    return added

# ---------------------------
# Multi-ligand helpers
# ---------------------------
def _find_ligand_sdfs(user_input: str | None,
                      glob_patterns=("ligand*.sdf", "ligands/*.sdf")) -> list[Path]:
    """
    Return a sorted list of SDFs to load.
    Priority:
      1) --input-ligand if provided (single-ligand mode)
      2) Any files matching glob patterns (multi-ligand mode)
      3) Fallback to ./ligand.sdf if present
    """
    if user_input:
        p = Path(user_input)
        return [p] if p.exists() else []
    found: list[Path] = []
    for pat in glob_patterns:
        for m in sorted(Path(".").glob(pat)):
            if m.is_file():
                found.append(m)
    if found:
        return found
    # fallback
    if Path("ligand.sdf").exists():
        return [Path("ligand.sdf")]
    return []

def _add_ligands(modeller, forcefield, sdf_paths: list[Path]) -> int:
    """
    Add one or more ligands from SDF paths to the modeller.
    Assign SMIRNOFF templates once for all of them.
    Rename residues to L01, L02, ... (3-char, PDB-safe).
    """
    if not sdf_paths:
        return 0
    off_mols: list[Molecule] = []
    for sdf in sdf_paths:
        lig = Molecule.from_file(str(sdf))
        if lig.n_conformers == 0:
            raise ValueError(f"Ligand has no conformer: {sdf}")
        lig.name = sdf.stem[:15]  # readable
        off_mols.append(lig)

    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    before_res = list(modeller.topology.residues())
    for lig in off_mols:
        lig_top = OFFTopology.from_molecules([lig]).to_openmm()
        lig_pos = to_openmm(lig.conformers[0])
        modeller.add(lig_top, lig_pos)

    # Rename newly added residues L01, L02, ...
    after_res = list(modeller.topology.residues())
    new_res = after_res[len(before_res):]
    for i, r in enumerate(new_res, start=1):
        try:
            r.name = f"L{i:02d}"  # PDB-safe 3-chars
        except Exception:
            pass
    return len(new_res)

# ---------------------------
# PDBx/mmCIF writer (for big systems)
# ---------------------------
def _write_cif(top, positions, out_path: Path):
    """
    Write PDBx/mmCIF to avoid PDB serial overflow (ChimeraX-friendly).
    """
    try:
        PDBxFile.writeFile(top, positions, str(out_path))
        print(f"🧾 Wrote mmCIF topology → {out_path}")
    except Exception as e:
        print(f"⚠️ Failed to write mmCIF ({out_path}): {e}")

# ---------------------------
# CLI
# ---------------------------
parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand(s), MixMD-from-Packmol).")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")

# Single or multi-ligand
parser.add_argument("--input-ligand", type=str, default=None,
                    help="Single ligand SDF (if omitted, auto-detect ligand*.sdf / ligands/*.sdf for multi-ligand MD)")

parser.add_argument("--n-steps", type=int, default=1_000_000, help="Number of MD steps (default: 1,000,000)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# MixMD-from-Packmol options
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0,
                    help="MD cubic box edge length; probes outside this cube are dropped.")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACEA,PHOL,ACOH",
                    help="Comma-separated probe residue names to import.")
parser.add_argument("--mixmd-from-packmol", action="store_true",
    help="Read placements CSV (build/<receptor>_mixmd_placements.csv) and place probes before solvation (implies --no-ligand).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
    help="Override path to placements CSV (default: build/<receptor>_mixmd_placements.csv).")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15,
    help="Safety margin from box edge when filtering centroids.")

# Hotspot logging
parser.add_argument("--hotspot-csv", type=str, default=None,
    help="If set with --mixmd-from-packmol, write probe centroids every --hotspot-stride to this CSV (default: mixmd_hotspots_no_ligand.csv).")
parser.add_argument("--hotspot-stride", type=int, default=1000,
    help="Stride (steps) for hotspot centroid logging.")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo (no ligands)
if getattr(args, "mixmd_from_packmol", False):
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# ---------------------------
# Load receptor
# ---------------------------
receptor_pdb = PDBFile(args.input_receptor)

# ---------------------------
# Force fields
# ---------------------------
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# ---------------------------
# Build modeller (ligand(s) or apo)
# ---------------------------
modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)

lig_added = 0
if not args.no_ligand and not args.mixmd_from_packmol:
    ligand_sdfs = _find_ligand_sdfs(args.input_ligand)
    if ligand_sdfs:
        print(f"🔗 Detected {len(ligand_sdfs)} ligand file(s): {[p.name for p in ligand_sdfs]}")
        lig_added = _add_ligands(modeller, forcefield, ligand_sdfs)
        print(f"✅ Added {lig_added} ligand residue(s).")
    else:
        print("ℹ️  No ligand SDF detected; proceeding with receptor-only (apo).")
else:
    if args.mixmd_from_packmol:
        print("ℹ️  --mixmd-from-packmol specified → running APO (ignoring ligands).")

# ---------------------------
# If MixMD requested, place probes before solvation
# ---------------------------
if args.mixmd_from_packmol:
    resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    added = _add_probes_from_packmol(
        modeller=modeller,
        forcefield=forcefield,
        receptor_path=receptor_path,
        sim_box_nm=float(args.mixmd_box_size_nm),
        edge_margin_nm=float(args.mixmd_edge_margin_nm),
        placements_csv=(Path(args.mixmd_placements_csv) if args.mixmd_placements_csv else None),
        resname_list=resnames
    )
    print(f"🧪 MixMD: added {added} probe instances.")

# ---------------------------
# Hydrogens
# ---------------------------
modeller.addHydrogens(forcefield)

with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

# ---------------------------
# Solvate
# ---------------------------
box_edge = float(args.mixmd_box_size_n m if args.mixmd_from_packmol else 7.0)
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(box_edge, box_edge, box_edge) * nanometer,
    ionicStrength=0.15 * molar,
    neutralize=True
)
# Write both PDB and mmCIF (CIF avoids PDB serial overflow in big boxes)
with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
_write_cif(modeller.topology, modeller.positions, Path(f"solvated_receptor_ligand{suffix}.cif"))
print(f"✅ Solvated system ready. Box = {box_edge:.2f} nm")

# ---------------------------
# Create System
# ---------------------------
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=HBonds
)

seed = int(args.seed)
barostat = MonteCarloBarostat(1 * bar, 300 * kelvin, 25)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
integrator.setRandomNumberSeed(seed)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# ---------------------------
# Hotspot reporter (MixMD only)
# ---------------------------
if args.mixmd_from_packmol:
    probe_res_hint = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    hotspot_csv = args.hotspot_csv or f"mixmd_hotspots{suffix}.csv"

    hrep = ProbeHotspotReporter(
        hotspot_csv,
        reportInterval=int(args.hotspot_stride),
        topology=simulation.topology,
        probe_resnames=probe_res_hint,
        center_mode="protein-com"
    )
    # Force an initial snapshot at step 0 so the CSV is never empty
    state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    hrep.report(simulation, state0)
    # Add to reporters BEFORE any steps so it logs during equilibration + production
    simulation.reporters.append(hrep)

# ---------------------------
# Run
# ---------------------------
print("🔹 Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter(f"minimized{suffix}.pdb", 100))

print("🔹 NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * kelvin, seed)
simulation.reporters.append(PDBReporter(f"nvt_equilibrated{suffix}.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(500)  # 1 ps

print("🔹 NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter(f"npt_equilibrated{suffix}.pdb", 500))
simulation.step(2500)  # 5 ps

print(f"🔥 Production MD: {args.n_steps:,} steps")
simulation.reporters.append(DCDReporter(f"production_md{suffix}.dcd", 1000))
simulation.reporters.append(StateDataReporter(f"production_md{suffix}.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(args.n_steps)

# Final snapshots (PDB + mmCIF)
with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
_write_cif(simulation.topology,
           simulation.context.getState(getPositions=True).getPositions(),
           Path(f"final_structure{suffix}.cif"))
print(f"🎉 MD complete → final_structure{suffix}.pdb / .cif saved.")
