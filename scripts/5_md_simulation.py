# 5_md_simulation.py
# Author: Iori Mochizuki + collaborator patch
# Updated: 2025-10-23
# Description: OpenMM MD (protein/RNA/ligand[s]) with optional MixMD-from-Packmol import.
# - Multi-ligand ready: auto-detects ligand*.sdf and/or ligands/*.sdf and adds all.
# - MixMD: reads Packmol placements CSV, adds probes pre-solvation, logs probe centroids each stride.
# - Hotspot CSV is receptor-centered (protein COM) and actually records beyond step 0 (robust schedule).
# - Writes both PDB and mmCIF (PDBx) to avoid 5-digit atom-serial issues in large boxes.

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

# =========================
# MixMD probe dictionary
# =========================
PROBE_LIBRARY = {
    # canonical 3-letter → SMILES
    "IPA":  "CC(C)O",           # isopropanol
    "ACN":  "CC#N",             # acetonitrile
    "IMD":  "c1[nH]cnc1",       # imidazole (neutral, N1–H)
    "ACM":  "CC(=O)N",          # acetamide
    "PHO":  "c1ccc(cc1)O",      # phenol
    "HAC":  "CC(=O)O",          # acetic acid
    # aliases accepted from older docs/runs
    "ACEA": "CC(=O)N",          # acetamide alias
    "PHOL": "c1ccc(cc1)O",      # phenol alias
    "ACOH": "CC(=O)O",          # acetic acid alias
}
PROBE_NAME_SET = set(PROBE_LIBRARY.keys())

def _default_packmol_csv(receptor_path: Path) -> Path:
    root = f"{receptor_path.stem}_mixmd"
    return Path("build") / f"{root}_placements.csv"


# =========================
# Hotspot Reporter (fixed)
# =========================
class ProbeHotspotReporter:
    """
    Write probe residue centroids (nm) every 'reportInterval' steps to CSV.
    Columns: step,resname,resid,x_nm,y_nm,z_nm

    Coordinates are stabilized in a receptor-centered frame:
       r' = r - COM_protein(t) + COM_protein(0)
    with COM based on protein CA atoms (fallback: protein heavy atoms).
    """

    def __init__(self, file_path, reportInterval, topology, probe_resnames=None, center_mode="protein-com"):
        self.file = open(file_path, "w")
        self.reportInterval = max(1, int(reportInterval))
        self.center_mode = center_mode
        self._have_ref_com = False
        self._ref_com = (0.0, 0.0, 0.0)

        # discover which probe residue names exist in topology
        present = {(res.name or "").upper().strip() for res in topology.residues()}
        auto_names = present & PROBE_NAME_SET
        if probe_resnames:
            user_names = {s.upper().strip() for s in probe_resnames}
            self.probe_res = auto_names | (user_names & PROBE_NAME_SET)
        else:
            self.probe_res = auto_names

        # index groups for probe residues + anchors for protein COM
        aa = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        groups, resid_map = [], []
        self._protein_anchor_indices = []

        for res in topology.residues():
            rn = (res.name or "").upper().strip()

            if rn in self.probe_res:
                idxs = [a.index for a in res.atoms()]
                if idxs:
                    groups.append(idxs)
                    resid_map.append((rn, int(res.id) if res.id is not None else 0))

            if rn in aa:
                # prefer CA atoms
                for a in res.atoms():
                    if (a.name or "").upper() == "CA":
                        self._protein_anchor_indices.append(a.index)

        if not self._protein_anchor_indices:
            # fallback: non-H protein atoms
            for res in topology.residues():
                rn = (res.name or "").upper().strip()
                if rn in aa:
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

    # *** Critical fix: robust schedule aligned to stride multiples ***
    def describeNextReport(self, simulation):
        # compute current integer step
        step = 0
        try:
            step = int(simulation.currentStep)
        except Exception:
            try:
                state = simulation.context.getState(getTime=True)
                t  = state.getTime()                       # has units
                dt = simulation.integrator.getStepSize()   # has units
                step = int(round((t/dt)))
            except Exception:
                step = 0
        n = self.reportInterval
        # number of steps until next multiple of n (never return 0)
        rem = step % n
        due = n if rem == 0 else (n - rem)
        return (due, True, False, False, False, True)

    @staticmethod
    def _com_of_indices(positions, idxs):
        sx = sy = sz = 0.0
        n = 0
        for i in idxs:
            v = positions[i]  # Vec3 (nm)
            sx += float(v.x); sy += float(v.y); sz += float(v.z)
            n += 1
        if n == 0:
            return (0.0, 0.0, 0.0)
        return (sx/n, sy/n, sz/n)

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=False)  # list of Vec3 (nm)

        # integer step again (for robust logging)
        try:
            step = int(simulation.currentStep)
        except Exception:
            try:
                t  = state.getTime()
                dt = simulation.integrator.getStepSize()
                step = int(round((t/dt)))
            except Exception:
                step = -1

        # receptor-centered recentering
        dx = dy = dz = 0.0
        if self.center_mode == "protein-com" and self._protein_anchor_indices:
            com_now = self._com_of_indices(pos, self._protein_anchor_indices)
            if not self._have_ref_com:
                self._ref_com = com_now
                self._have_ref_com = True
            dx = self._ref_com[0] - com_now[0]
            dy = self._ref_com[1] - com_now[1]
            dz = self._ref_com[2] - com_now[2]

        # write all probe centroids
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


# =========================
# Helpers: probes from Packmol
# =========================
def _openff_conf_to_angstrom_list(off_mol: Molecule):
    """Return [(x,y,z) in Å] from an OpenFF Molecule conformer (no NumPy)."""
    conf = off_mol.conformers[0]  # OpenFF quantity with units
    try:
        arr = conf.to("angstrom").m  # N×3 numeric via Pint
        return [(float(r[0]), float(r[1]), float(r[2])) for r in arr]
    except Exception:
        q = to_openmm(conf)  # Quantity(list(Vec3), angstrom)
        out = []
        for v in q.value_in_unit(angstrom):
            out.append((float(v[0]), float(v[1]), float(v[2])))
        return out

def _positions_for_centroid_nm(off_mol: Molecule, centroid_nm_tuple):
    """Return [Vec3]*nm for `off_mol`, centered at centroid_nm_tuple (x,y,z in nm)."""
    coords_A = _openff_conf_to_angstrom_list(off_mol)
    n = len(coords_A)
    if n == 0:
        raise ValueError("OFF molecule has no coordinates.")
    sx = sy = sz = 0.0
    for xA, yA, zA in coords_A:
        sx += xA; sy += yA; sz += zA
    cx = sx/n; cy = sy/n; cz = sz/n
    tx, ty, tz = (float(centroid_nm_tuple[0]), float(centroid_nm_tuple[1]), float(centroid_nm_tuple[2]))
    out = []
    for xA, yA, zA in coords_A:
        out.append(mm.Vec3((xA-cx)*0.1 + tx, (yA-cy)*0.1 + ty, (zA-cz)*0.1 + tz))
    return [v*nanometer for v in out]

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

def _add_probes_from_packmol(modeller: Modeller,
                             forcefield: ForceField,
                             receptor_path: Path,
                             sim_box_nm: float,
                             edge_margin_nm: float,
                             placements_csv: Path | None = None,
                             resname_list: list[str] | None = None) -> int:
    """
    Add probes at Packmol centroids cropped to the intended MD box.
    Renames newly added residues to the probe resname so the reporter can find them.
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
        m = Molecule.from_smiles(PROBE_LIBRARY[res], allow_undefined_stereo=True)
        m.generate_conformers(n_conformers=1)
        m.name = res
        off_mols.append(m)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    top_cache = {m.name: OFFTopology.from_molecules([m]).to_openmm() for m in off_mols}
    mol_cache = {m.name: m for m in off_mols}

    def _res_list(): return list(modeller.topology.residues())

    added = 0
    for (resname, pack_resid, cen_nm) in rows:
        top = top_cache[resname]
        off = mol_cache[resname]
        pos = _positions_for_centroid_nm(off, cen_nm)

        before = len(_res_list())
        modeller.add(top, pos)
        new_res = _res_list()[before:]

        for r in new_res:
            try: r.name = resname
            except Exception: pass
            try: r.id = str(pack_resid)
            except Exception: pass

        added += 1

    print(f"✅ Placed {added} probe molecules (cropped to {sim_box_nm:.1f} nm box).")
    return added


# =========================
# Helpers: multi‑ligand load
# =========================
def _find_ligand_sdfs(patterns=("ligand*.sdf",), folders=("ligands", "Ligands", "LIGANDS")):
    found = set()
    for pat in patterns:
        for p in glob.glob(pat):
            if p.lower().endswith(".sdf"):
                found.add(Path(p).resolve())
    for folder in folders:
        d = Path(folder)
        if d.is_dir():
            for p in d.glob("*.sdf"):
                found.add(p.resolve())
    # deterministic order
    return sorted(found, key=lambda p: p.name)

def _add_ligands(modeller: Modeller, forcefield: ForceField, sdf_paths: list[Path]) -> int:
    """
    Add one or many ligands (already positioned in SDF coords).
    Each SDF should contain a single molecule with 3D conformer.
    """
    if not sdf_paths:
        return 0
    off_mols = []
    tops = []
    poss = []
    for sdf in sdf_paths:
        lig = Molecule.from_file(str(sdf))
        if len(lig.conformers) == 0:
            # attempt one conformer (keeps absolute coords if present)
            lig.generate_conformers(n_conformers=1)
        off_mols.append(lig)
        tops.append(OFFTopology.from_molecules([lig]).to_openmm())
        poss.append(to_openmm(lig.conformers[0]))
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    for top, pos in zip(tops, poss):
        modeller.add(top, pos)
    print(f"✅ Added {len(sdf_paths)} ligand(s): " + ", ".join(p.name for p in sdf_paths))
    return len(sdf_paths)


# =========================
# CLI
# =========================
parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand[s], MixMD-from-Packmol).")

# receptor / ligand
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default=None,
                    help="Single-ligand SDF (kept for backward compatibility).")
parser.add_argument("--ligand-pattern", type=str, default="ligand*.sdf",
                    help="Glob pattern to pick up multiple SDF ligands (default: ligand*.sdf).")
parser.add_argument("--ligand-dir", type=str, default=None,
                    help="Optional folder to scan for *.sdf (e.g., ligands/).")

# MixMD
parser.add_argument("--mixmd-from-packmol", action="store_true",
    help="Read Packmol placements CSV (build/<receptor>_mixmd_placements.csv) and place probes before solvation (implies --no-ligand).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
    help="Override path to placements CSV (default: build/<receptor>_mixmd_placements.csv).")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACM,PHO,HAC,ACEA,PHOL,ACOH",
    help="Comma-separated probe residue names to import.")
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0,
    help="MD cubic box edge length; probes outside this cube are dropped.")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15,
    help="Safety margin from box edge when filtering centroids.")

# MD
parser.add_argument("--n-steps", type=int, default=1_000_000, help="Production MD steps (default: 1,000,000)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# Hotspots
parser.add_argument("--hotspot-csv", type=str, default=None,
    help="If set with --mixmd-from-packmol, write probe centroids every --hotspot-stride to this CSV (default: mixmd_hotspots_no_ligand.csv).")
parser.add_argument("--hotspot-stride", type=int, default=1000,
    help="Stride (steps) for hotspot centroid logging.")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo
if getattr(args, "mixmd_from_packmol", False):
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# =========================
# Load receptor & forcefield
# =========================
receptor_pdb = PDBFile(args.input_receptor)

if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# =========================
# Build modeller (ligands or apo)
# =========================
modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)

added_ligands = 0
if not args.no_ligand:
    # gather SDFs
    candidates = []
    if args.input_ligand:
        p = Path(args.input_ligand)
        if p.exists() and p.suffix.lower()==".sdf":
            candidates.append(p.resolve())
    # glob pattern in CWD
    for p in _find_ligand_sdfs(patterns=(args.ligand_pattern,)):
        if p not in candidates:
            candidates.append(p)
    # optional folder
    if args.ligand_dir:
        d = Path(args.ligand_dir)
        if d.is_dir():
            for p in d.glob("*.sdf"):
                rp = p.resolve()
                if rp not in candidates:
                    candidates.append(rp)

    added_ligands = _add_ligands(modeller, forcefield, candidates)
    if added_ligands == 0 and args.input_ligand:
        # Back-compat: try single path even if extension is different
        lig = Molecule.from_file(args.input_ligand)
        if len(lig.conformers) == 0:
            lig.generate_conformers(n_conformers=1)
        smirnoff = SMIRNOFFTemplateGenerator(molecules=[lig])
        forcefield.registerTemplateGenerator(smirnoff.generator)
        modeller.add(OFFTopology.from_molecules([lig]).to_openmm(), to_openmm(lig.conformers[0]))
        added_ligands = 1
    print(f"✅ Ligand mode: added {added_ligands} ligand(s).")
else:
    print("✅ No ligand: receptor-only (apo) setup")

# =========================
# MixMD probes (before solvation)
# =========================
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

# Hydrogens after all non-water species present
modeller.addHydrogens(forcefield)

with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else f'holo: {added_ligands} ligand(s)'}).")

# =========================
# Solvate (box size)
# =========================
box_edge = float(args.mixmd_box_size_nm if args.mixmd_from_packmol else 7.0)
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(box_edge, box_edge, box_edge) * nanometer,
    ionicStrength=0.15 * molar,
    neutralize=True
)
with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ Solvated system ready. Box = {box_edge:.2f} nm")

# =========================
# System & Integrator
# =========================
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

# =========================
# Hotspot reporter
# =========================
if args.mixmd_from_packmol:
    probe_resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    hotspot_csv = args.hotspot_csv or f"mixmd_hotspots{suffix}.csv"

    hrep = ProbeHotspotReporter(
        hotspot_csv,
        reportInterval=int(args.hotspot_stride),
        topology=simulation.topology,
        probe_resnames=probe_resnames,
        center_mode="protein-com"
    )

    # Force an initial snapshot (step 0), then append reporter for scheduled logging
    state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    hrep.report(simulation, state0)
    simulation.reporters.append(hrep)

# =========================
# Run MD
# =========================
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
# --- DCD reference structure (write at the moment we start DCD) ---
# This snapshot has exactly the same atom order/count ChimeraX will expect.
prod_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
ref_pdb = f"production_start_structure{suffix}.pdb"
ref_cif = f"production_start_structure{suffix}.cif"

with open(ref_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, prod_state.getPositions(), f)

# Keep mmCIF too (useful when atom serials exceed 99,999)
with open(ref_cif, "w") as f:
    PDBxFile.writeFile(simulation.topology, prod_state.getPositions(), f)

print(f"📌 DCD reference written → {ref_pdb} (and {ref_cif})")
# ------------------------------------------------------------------
simulation.reporters.append(DCDReporter(f"production_md{suffix}.dcd", 1000))
simulation.reporters.append(StateDataReporter(f"production_md{suffix}.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(args.n_steps)

# extra: ensure hotspot reporter at final state too
if args.mixmd_from_packmol:
    try:
        final_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        for rep in simulation.reporters:
            if isinstance(rep, ProbeHotspotReporter):
                rep.report(simulation, final_state)
    except Exception:
        pass

# =========================
# Save final structures (PDB + mmCIF)
# =========================
final_positions = simulation.context.getState(getPositions=True).getPositions()
with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, final_positions, f)
print(f"✅ PDB written → final_structure{suffix}.pdb")

# mmCIF avoids 5-digit serial overflow; very important for large boxes
with open(f"final_structure{suffix}.cif", "w") as f:
    PDBxFile.writeFile(simulation.topology, final_positions, f)
print(f"✅ mmCIF written → final_structure{suffix}.cif")

print("🎉 MD complete.")
