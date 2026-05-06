#!/usr/bin/env python3
"""
pada1_md_pipeline.py
====================
Full MD pipeline for PaDa-I UPO (PDB: 5OXU) with explicit heme parameterization.

Pipeline stages:
  1. PDB preparation      — download 5OXU, clean, protonate, split
  2. MCPB.py step 1       — build QM model, generate Gaussian inputs
  3. Gaussian             — geometry opt + force constants + ESP charges
  4. MCPB.py steps 2–4   — extract FF params, generate tleap input
  5. tleap                — assemble solvated system → prmtop/inpcrd
  6. OpenMM via UnoMD     — energy minimisation → NVT equilibration → production MD
  7. Analysis             — RMSD, RMSF, Fe–S bond, Phe191 dihedral, channel width

Usage:
  python pada1_md_pipeline.py [--stage N] [--config config.yaml] [--ligand ligand.sdf]

Requirements:
  conda env:  ambertools  openmm  mdanalysis  pdbfixer  requests
  pip:        unomd
  binary:     g16 (Gaussian 16) in PATH
  optional:   cuda toolkit for GPU acceleration

Author: BioBuilders 2026
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import yaml

# ── optional heavy imports (only needed in later stages) ─────────────────────
def _import_openmm():
    from openmm.app import (AmberInpcrdFile, AmberPrmtopFile, DCDReporter,
                             HBonds, PME, Simulation, StateDataReporter)
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
    return locals()

def _import_mda():
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, dihedrals
    from MDAnalysis.analysis.distances import dist
    return mda, rms, dihedrals, dist

# ── logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("pada1_pipeline")


# ═════════════════════════════════════════════════════════════════════════════
#  DEFAULT CONFIGURATION
# ═════════════════════════════════════════════════════════════════════════════

DEFAULT_CONFIG = {
    # ── paths ────────────────────────────────────────────────────────────────
    "pdb_id":        "5OXU",
    "workdir":       "pada1_md",
    "ligand_sdf":    None,            # optional: path to substrate .sdf file

    # ── structure preparation ─────────────────────────────────────────────────
    "ph":            7.0,
    "keep_glycans":  False,           # False = strip GlcNAc for cleaner first run
    "mg_resname":    "MG",            # structural Mg2+ residue name in PDB
    "heme_resname":  "HEM",
    "cys_axial":     36,              # residue number of the Cys axial ligand

    # ── MCPB.py ──────────────────────────────────────────────────────────────
    "mcpb_cutoff":   2.8,             # metal coordination cutoff (Å)
    "mcpb_ff":       "ff14SB",
    "mcpb_gaff":     "gaff2",

    # ── QM engine ────────────────────────────────────────────────────────────
    "qm_engine":     "orca",          # "orca" (free, academic) or "gaussian"
    "orca_exe":      "/root/orca_6_1_1/orca",
    "orca_nproc":    4,               # parallel cores for ORCA %pal
    "orca_mem_mb":   4000,            # MB per core for ORCA %maxcore
    "gaussian_exe":  "g16",           # used only if qm_engine = "gaussian"
    "gaussian_nproc": 8,
    "gaussian_mem":  "16GB",
    "spin_mult":     6,               # ferric high-spin S=5/2 → multiplicity=6
    "charge":        0,               # overall QM model charge

    # ── tleap ────────────────────────────────────────────────────────────────
    "solv_buffer":   12.0,            # Å padding around protein
    "ion_conc":      0.15,            # physiological NaCl (mol/L)
    "water_model":   "TIP3PBOX",

    # ── MD settings ──────────────────────────────────────────────────────────
    "platform":      "CUDA",          # CUDA | CPU | OpenCL
    "precision":     "mixed",
    "temperature_K": 300,
    "friction":      1.0,             # ps^-1
    "timestep_fs":   2.0,
    "emin_tol":      5.0,             # kJ/mol/nm
    "equil_ns":      1.0,             # NVT equilibration length
    "prod_ns":       100.0,           # production MD length
    "save_interval": 10000,           # frames every N steps (~20 ps at 2 fs)
    "npt_prod":      True,            # NPT for production

    # ── analysis ─────────────────────────────────────────────────────────────
    "ref_pdb":       "5OXU",          # reference structure for RMSD
    "fe_s_warn_ang": 2.8,             # Å threshold for Fe–S bond drift warning
    "channel_res":   [76, 191],       # CA pair defining channel entrance width
}


# ═════════════════════════════════════════════════════════════════════════════
#  UTILITIES
# ═════════════════════════════════════════════════════════════════════════════

def load_config(path: str | None) -> dict:
    cfg = DEFAULT_CONFIG.copy()
    if path and Path(path).exists():
        with open(path) as f:
            overrides = yaml.safe_load(f)
        if overrides:
            cfg.update(overrides)
        log.info(f"Config loaded from {path}")
    return cfg


def run(cmd: str, cwd: str | None = None, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command, stream stdout, raise on failure."""
    log.info(f"$ {cmd}")
    result = subprocess.run(
        cmd, shell=True, cwd=cwd,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if result.stdout:
        for line in result.stdout.splitlines():
            log.debug(line)
    if check and result.returncode != 0:
        log.error(f"Command failed (exit {result.returncode}):\n{result.stdout[-2000:]}")
        raise RuntimeError(f"Command failed: {cmd}")
    return result


def workpath(cfg: dict, *parts) -> Path:
    p = Path(cfg["workdir"]).joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def require_binary(name: str):
    if not (shutil.which(name) or (Path(name).is_absolute() and Path(name).exists())):
        raise EnvironmentError(
            f"'{name}' not found. "
            f"Install AmberTools (conda install -c conda-forge ambertools), "
            f"or set the correct path in your config (orca_exe / gaussian_exe)."
        )


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 1 — PDB PREPARATION
# ═════════════════════════════════════════════════════════════════════════════

def stage1_prepare_pdb(cfg: dict, paths: dict = None) -> dict:
    """
    1a. Download 5OXU from RCSB.
    1b. Clean with pdb4amber: protonate, renumber, flag CYM.
    1c. Split into protein / heme / ligand files.
    Returns dict of output file paths.
    """
    log.info("=" * 60)
    log.info("STAGE 1 — PDB PREPARATION")
    log.info("=" * 60)

    require_binary("pdb4amber")

    wdir = Path(cfg["workdir"]) / "01_prep"
    wdir.mkdir(parents=True, exist_ok=True)

    pdb_id   = cfg["pdb_id"]
    raw_pdb  = wdir / f"{pdb_id}_raw.pdb"
    prep_pdb = wdir / f"{pdb_id}_prep.pdb"
    prot_pdb = wdir / f"{pdb_id}_protein.pdb"
    heme_pdb = wdir / f"{pdb_id}_heme.pdb"
    mg_pdb   = wdir / f"{pdb_id}_mg.pdb"

    # ── 1a: download ─────────────────────────────────────────────────────────
    if not raw_pdb.exists():
        log.info(f"Downloading {pdb_id} from RCSB...")
        import urllib.request
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, raw_pdb)
        log.info(f"  Saved to {raw_pdb}")
    else:
        log.info(f"  Using cached {raw_pdb}")

    # ── 1b: clean with pdb4amber ─────────────────────────────────────────────
    log.info("Running pdb4amber...")
    run(
        f"pdb4amber -i {raw_pdb.name} -o {prep_pdb.name} "
        f"--reduce --add-missing-atoms --no-conect",
        cwd=str(wdir.resolve())
    )

    # ── 1c: parse and split ───────────────────────────────────────────────────
    log.info("Splitting PDB into protein / heme / Mg...")
    lines = raw_pdb.read_text().splitlines()

    heme_resname = cfg["heme_resname"]
    mg_resname   = cfg["mg_resname"]
    cys_axial    = cfg["cys_axial"]
    keep_glycans = cfg["keep_glycans"]

    protein_lines, heme_lines, mg_lines = [], [], []

    # altloc filter: keep blank or 'A' only
    altloc_re = re.compile(r"^(ATOM|HETATM).{10}([A-Z ])")

    for line in lines:
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        m = altloc_re.match(line)
        if m and m.group(2) not in (" ", "A"):
            continue  # skip altloc B/C

        resname = line[17:20].strip()
        resnum  = int(line[22:26].strip())

        if rec == "HETATM" and resname == heme_resname:
            heme_lines.append(line)
        elif rec == "HETATM" and resname == mg_resname:
            mg_lines.append(line)
        elif rec == "HETATM" and resname in ("NAG", "FUC", "MAN") and not keep_glycans:
            pass  # strip glycans
        elif rec in ("ATOM", "HETATM"):
            # rename Cys36 → CYM (deprotonated thiolate axial ligand)
            if rec == "ATOM" and resname == "CYS" and resnum == cys_axial:
                line = line[:17] + "CYM" + line[20:]
            protein_lines.append(line)

    prot_pdb.write_text("\n".join(protein_lines) + "\nEND\n")
    heme_pdb.write_text("\n".join(heme_lines) + "\nEND\n")
    mg_pdb.write_text("\n".join(mg_lines) + "\nEND\n")

    log.info(f"  Protein: {len(protein_lines)} atoms → {prot_pdb}")
    log.info(f"  Heme:    {len(heme_lines)} atoms → {heme_pdb}")
    log.info(f"  Mg:      {len(mg_lines)} atoms → {mg_pdb}")

    # ── 1d: validate Cys36 rename ────────────────────────────────────────────
    cym_count = sum(1 for l in protein_lines if "CYM" in l[17:20])
    if cym_count == 0:
        raise ValueError(
            f"Cys{cys_axial} was not renamed to CYM. "
            f"Check cys_axial in config (currently {cys_axial})."
        )
    log.info(f"  CYM axial ligand confirmed at residue {cys_axial} ({cym_count} atoms).")

    return {
        "prep_pdb":   str(prep_pdb),
        "prot_pdb":   str(prot_pdb),
        "heme_pdb":   str(heme_pdb),
        "mg_pdb":     str(mg_pdb),
        "wdir":       str(wdir),
    }


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 2 — MCPB.py STEP 1 (build QM model)
# ═════════════════════════════════════════════════════════════════════════════

def stage2_mcpb_step1(cfg: dict, paths: dict) -> dict:
    """
    Merge protein + heme PDB, write MCPB.py input file, run step 1.
    Outputs Gaussian .com files.
    """
    log.info("=" * 60)
    log.info("STAGE 2 — MCPB.py STEP 1")
    log.info("=" * 60)

    require_binary("MCPB.py")

    wdir = Path(cfg["workdir"]) / "02_mcpb"
    wdir.mkdir(parents=True, exist_ok=True)

    pdb_id    = cfg["pdb_id"].lower()
    merged    = wdir / f"{pdb_id}_merged.pdb"
    mcpb_in   = wdir / f"{pdb_id}_heme.in"

    # ── merge protein + heme ─────────────────────────────────────────────────
    log.info("Merging protein + heme PDB...")
    with open(merged, "w") as out:
        for src in [paths["prot_pdb"], paths["heme_pdb"]]:
            for line in Path(src).read_text().splitlines():
                if line.strip() not in ("END", ""):
                    out.write(line + "\n")
    out_text = merged.read_text()
    # find the Fe atom serial number
    fe_serial = None
    for line in out_text.splitlines():
        if line[:6].strip() == "HETATM" and line[12:14].strip() == "FE":
            fe_serial = int(line[6:11].strip())
            break
    if fe_serial is None:
        raise ValueError("Fe atom not found in merged PDB. Check heme PDB content.")
    log.info(f"  Fe atom serial number: {fe_serial}")

    # ── write MCPB.py input ──────────────────────────────────────────────────
    qm_engine = cfg.get("qm_engine", "orca").lower()
    orca_block = (
        f"\nquantum_soft orca\norca_mem_mb {cfg.get('orca_mem_mb', 4000)}\n"
        if qm_engine == "orca" else ""
    )
    log.info(f"Writing MCPB.py input → {mcpb_in}")
    mcpb_in.write_text(textwrap.dedent(f"""
        original_pdb {merged.name}
        group_name {pdb_id}_heme
        cut_off {cfg['mcpb_cutoff']}
        ion_ids {fe_serial}
        ion_mol2files FE.mol2
        naa_mol2files HEM.mol2
        frcmod_files HEM.frcmod
        ff_choice {cfg['mcpb_ff']}
        gaff_version {cfg['mcpb_gaff']}
        large_opt 1
        conf_search 0
    """).strip() + orca_block + "\n")

    # ── write a minimal Fe mol2 ───────────────────────────────────────────────
    fe_mol2 = wdir / "FE.mol2"
    if not fe_mol2.exists():
        fe_mol2.write_text(textwrap.dedent("""
            @<TRIPOS>MOLECULE
            FE
             1  0  0  0  0
            SMALL
            NO_CHARGES

            @<TRIPOS>ATOM
                  1 FE          0.0000    0.0000    0.0000 Fe      1  FE        3.0000
            @<TRIPOS>BOND
        """).strip() + "\n")
        log.info(f"  Wrote minimal Fe.mol2 → {fe_mol2}")

    # ── check for pre-built Shahrokh heme params ──────────────────────────────
    hem_frcmod = wdir / "HEM.frcmod"
    hem_mol2   = wdir / "HEM.mol2"
    if not hem_frcmod.exists() or not hem_mol2.exists():
        log.warning(
            "HEM.frcmod and/or HEM.mol2 not found in 02_mcpb/.\n"
            "  → Download Shahrokh et al. 2012 parameters from:\n"
            "    https://pmc.ncbi.nlm.nih.gov/articles/PMC3242737/\n"
            "    (Supplementary data: NIHMS316516-supplement.doc)\n"
            "  → Place HEM.frcmod and HEM.mol2 in: " + str(wdir)
        )

    # ── run MCPB.py step 1 ───────────────────────────────────────────────────
    # MCPB.py generates .inp for ORCA (quantum_soft orca), .com for Gaussian
    ext = "inp" if qm_engine == "orca" else "com"
    log.info(f"Running MCPB.py step 1 (QM engine: {qm_engine})...")
    run(f"MCPB.py -i {mcpb_in.name} -s 1", cwd=str(wdir))

    small_opt = wdir / f"{pdb_id}_heme_small_opt.{ext}"
    small_fc  = wdir / f"{pdb_id}_heme_small_fc.{ext}"
    large_mk  = wdir / f"{pdb_id}_heme_large_mk.{ext}"

    for f in [small_opt, small_fc, large_mk]:
        if not f.exists():
            raise FileNotFoundError(f"Expected QM input not generated: {f}")
        log.info(f"  Generated: {f.name}")

    spin   = cfg["spin_mult"]
    charge = cfg["charge"]

    if qm_engine == "gaussian":
        nproc = cfg["gaussian_nproc"]
        mem   = cfg["gaussian_mem"]
        for com_file in [small_opt, small_fc, large_mk]:
            txt = com_file.read_text()
            txt = re.sub(r"%[Nn]proc\w*=\d+", f"%nproc={nproc}", txt)
            txt = re.sub(r"%[Mm]em=\S+",      f"%mem={mem}",     txt)
            txt = re.sub(r"^(\s*0\s+1\s*)$", f"  {charge} {spin}", txt, flags=re.MULTILINE)
            com_file.write_text(txt)
        log.info(f"  Patched Gaussian inputs: charge={charge}, mult={spin}, nproc={nproc}, mem={mem}")
    # ORCA: MCPB.py writes correct ORCA format directly — no patching needed

    return {
        "wdir":       str(wdir),
        "merged_pdb": str(merged),
        "mcpb_in":    str(mcpb_in),
        "small_opt":  str(small_opt),
        "small_fc":   str(small_fc),
        "large_mk":   str(large_mk),
        "pdb_id":     pdb_id,
        "qm_engine":  qm_engine,
    }


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 3 — QM JOBS (ORCA or Gaussian)
# ═════════════════════════════════════════════════════════════════════════════

def stage3_gaussian(cfg: dict, paths: dict) -> dict:
    """
    Submit three QM jobs in sequence using ORCA or Gaussian.
    Jobs must complete before MCPB.py step 2 can proceed.
    """
    qm_engine = paths.get("qm_engine", cfg.get("qm_engine", "gaussian")).lower()

    log.info("=" * 60)
    log.info(f"STAGE 3 — QM CALCULATIONS ({qm_engine.upper()})")
    log.info("=" * 60)

    wdir   = Path(paths["wdir"])
    pdb_id = paths["pdb_id"]

    if qm_engine == "orca":
        return _stage3_orca(cfg, paths, wdir, pdb_id)
    else:
        return _stage3_gaussian(cfg, paths, wdir, pdb_id)


def _stage3_gaussian(cfg, paths, wdir, pdb_id):
    exe = cfg["gaussian_exe"]
    require_binary(exe)

    jobs = [
        (paths["small_opt"], f"{pdb_id}_heme_small_opt.log",
         "Geometry optimisation of Fe coordination shell"),
        (paths["small_fc"],  f"{pdb_id}_heme_small_fc.log",
         "Force constant scan (Fe–S, Fe–N bonds/angles)"),
        (paths["large_mk"],  f"{pdb_id}_heme_large_mk.log",
         "Full-model MK ESP charges (RESP fitting)"),
    ]

    log_files = {}
    for com, logname, desc in jobs:
        logpath = wdir / logname
        if logpath.exists() and _qm_converged_gaussian(logpath):
            log.info(f"  Skipping {logname} — already converged.")
        else:
            log.info(f"  Running: {desc}")
            run(f"{exe} < {Path(com).name} > {logname}", cwd=str(wdir))
            if not _qm_converged_gaussian(logpath):
                raise RuntimeError(
                    f"Gaussian job did not converge: {logpath}\n"
                    f"Check {logname} for SCF convergence / spin state errors."
                )
            log.info(f"  Converged: {logname}")
        log_files[Path(com).stem] = str(logpath)

    return {**paths, "qm_logs": log_files}


def _stage3_orca(cfg, paths, wdir, pdb_id):
    exe = cfg.get("orca_exe", "orca")
    require_binary(exe)

    jobs = [
        (paths["small_opt"], f"{pdb_id}_heme_small_opt.out",
         "Geometry optimisation of Fe coordination shell"),
        (paths["small_fc"],  f"{pdb_id}_heme_small_fc.out",
         "Force constant scan (Fe–S, Fe–N bonds/angles)"),
        (paths["large_mk"],  f"{pdb_id}_heme_large_mk.out",
         "Full-model ESP charges (RESP fitting)"),
    ]

    log_files = {}
    for inp, outname, desc in jobs:
        outpath = wdir / outname
        if outpath.exists() and _qm_converged_orca(outpath):
            log.info(f"  Skipping {outname} — already converged.")
        else:
            log.info(f"  Running: {desc}")
            run(f"{exe} {Path(inp).name} > {outname}", cwd=str(wdir))
            if not _qm_converged_orca(outpath):
                raise RuntimeError(
                    f"ORCA job did not converge: {outpath}\n"
                    f"Check {outname} for SCF convergence / spin state errors."
                )
            log.info(f"  Converged: {outname}")
        log_files[Path(inp).stem] = str(outpath)

    return {**paths, "qm_logs": log_files}


def _qm_converged_gaussian(logpath: Path) -> bool:
    if not logpath.exists():
        return False
    return "Normal termination" in logpath.read_text()[-3000:]


def _qm_converged_orca(logpath: Path) -> bool:
    if not logpath.exists():
        return False
    return "ORCA TERMINATED NORMALLY" in logpath.read_text()[-3000:]


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 4 — MCPB.py STEPS 2–4 (extract FF params + tleap input)
# ═════════════════════════════════════════════════════════════════════════════

def stage4_mcpb_steps234(cfg: dict, paths: dict) -> dict:
    log.info("=" * 60)
    log.info("STAGE 4 — MCPB.py STEPS 2–4")
    log.info("=" * 60)

    require_binary("MCPB.py")

    wdir   = Path(paths["wdir"])
    mcpb_in = Path(paths["mcpb_in"]).name
    pdb_id  = paths["pdb_id"]

    for step_num, desc in [
        (2, "Extract bonded parameters from Gaussian force constant scan"),
        (3, "Fit RESP charges from MK ESP log"),
        (4, "Generate tleap input script"),
    ]:
        log.info(f"  MCPB.py step {step_num}: {desc}")
        run(f"MCPB.py -i {mcpb_in} -s {step_num}", cwd=str(wdir))

    frcmod_out = wdir / f"{pdb_id}_heme_small.frcmod"
    resp_mol2  = wdir / f"{pdb_id}_heme_large_RESP.mol2"
    tleap_in   = wdir / f"{pdb_id}_heme_tleap.in"

    for f in [frcmod_out, resp_mol2, tleap_in]:
        if not f.exists():
            raise FileNotFoundError(f"MCPB.py did not produce expected output: {f}")
        log.info(f"  Generated: {f.name}")

    return {**paths, "frcmod": str(frcmod_out), "resp_mol2": str(resp_mol2), "tleap_in": str(tleap_in)}


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 5 — tleap SYSTEM ASSEMBLY
# ═════════════════════════════════════════════════════════════════════════════

def stage5_tleap(cfg: dict, paths: dict) -> dict:
    log.info("=" * 60)
    log.info("STAGE 5 — TLEAP SYSTEM ASSEMBLY")
    log.info("=" * 60)

    require_binary("tleap")

    wdir   = Path(cfg["workdir"]) / "03_tleap"
    wdir.mkdir(parents=True, exist_ok=True)

    pdb_id    = paths["pdb_id"]
    prot_pdb  = Path(paths["prot_pdb"]).resolve()
    frcmod    = Path(paths["frcmod"]).resolve()
    resp_mol2 = Path(paths["resp_mol2"]).resolve()
    hem_frcmod = (Path(paths["wdir"]) / "HEM.frcmod").resolve()

    prmtop  = wdir / f"{pdb_id}_solvated.prmtop"
    inpcrd  = wdir / f"{pdb_id}_solvated.inpcrd"
    dry_pdb = wdir / f"{pdb_id}_dry.pdb"

    # ── optional ligand parameterization ─────────────────────────────────────
    lig_block = ""
    if cfg.get("ligand_sdf"):
        lig_frcmod, lig_mol2 = _parameterize_ligand(cfg, wdir)
        lig_block = textwrap.dedent(f"""
            # Substrate ligand
            loadamberparams {lig_frcmod}
            LIG = loadmol2 {lig_mol2}
        """)

    # ── write tleap input ────────────────────────────────────────────────────
    cys_axial = cfg["cys_axial"]
    buf       = cfg["solv_buffer"]
    water     = cfg["water_model"]
    mg_res    = cfg.get("mg_resname", "MG")

    tleap_script = textwrap.dedent(f"""
        source leaprc.protein.ff14SB
        source leaprc.water.tip3p
        source leaprc.gaff2

        # MCPB-derived heme bonded parameters (Fe-S + Fe-N force constants)
        loadamberparams {frcmod}

        # Shahrokh 2012 heme macrocycle parameters (porphyrin ring, propionates)
        loadamberparams {hem_frcmod}

        # Full heme + CYM coordination shell with RESP charges
        HEM = loadmol2 {resp_mol2}

        {lig_block}

        # Load protein (Cys{cys_axial} = CYM, deprotonated thiolate)
        mol = loadpdb {prot_pdb}

        # Explicit Fe-S and Fe-N4 bonds (bonded model for heme iron)
        bond mol.{cys_axial}.SG  mol.HEM.FE
        bond mol.HEM.NA          mol.HEM.FE
        bond mol.HEM.NB          mol.HEM.FE
        bond mol.HEM.NC          mol.HEM.FE
        bond mol.HEM.ND          mol.HEM.FE

        # Check structure before solvation
        check mol

        # Neutralise then add physiological NaCl
        addions mol Na+ 0
        addions mol Cl- 0

        # Solvate in TIP3P water box
        solvatebox mol {water} {buf}

        # Save dry structure for reference
        savepdb mol {dry_pdb}

        # Save AMBER topology and coordinates
        saveamberparm mol {prmtop} {inpcrd}

        quit
    """).strip()

    tleap_file = wdir / "tleap.in"
    tleap_file.write_text(tleap_script)
    log.info(f"  tleap script → {tleap_file}")

    run(f"tleap -f {tleap_file}", cwd=str(wdir))

    for f in [prmtop, inpcrd]:
        if not f.exists():
            raise FileNotFoundError(f"tleap did not produce: {f}")
        log.info(f"  Generated: {f.name}  ({f.stat().st_size / 1e6:.1f} MB)")

    # ── quick sanity check: count Fe bonds in prmtop ──────────────────────────
    prmtop_txt = prmtop.read_text(errors="replace")
    fe_bonds = prmtop_txt.count("FE")
    log.info(f"  Fe occurrences in prmtop: {fe_bonds}  (expect ≥5 for 1×Fe–S + 4×Fe–N)")
    if fe_bonds < 5:
        log.warning("  Fe bond count looks low — check tleap output for bond errors.")

    return {**paths, "prmtop": str(prmtop), "inpcrd": str(inpcrd), "dry_pdb": str(dry_pdb)}


def _parameterize_ligand(cfg: dict, wdir: Path) -> tuple[str, str]:
    """Run antechamber + parmchk2 to parameterize the substrate ligand."""
    require_binary("antechamber")
    require_binary("parmchk2")

    sdf = cfg["ligand_sdf"]
    lig_mol2   = wdir / "ligand.mol2"
    lig_frcmod = wdir / "ligand.frcmod"

    if not lig_mol2.exists():
        log.info(f"  Parameterizing ligand from {sdf}...")
        run(
            f"antechamber -i {sdf} -fi sdf "
            f"-o {lig_mol2} -fo mol2 "
            f"-c bcc -s 2 -nc 0",
            cwd=str(wdir)
        )
    if not lig_frcmod.exists():
        run(f"parmchk2 -i {lig_mol2} -f mol2 -o {lig_frcmod}", cwd=str(wdir))

    log.info(f"  Ligand params: {lig_mol2.name}, {lig_frcmod.name}")
    return str(lig_frcmod), str(lig_mol2)


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 6 — OpenMM MD (via UnoMD orchestration)
# ═════════════════════════════════════════════════════════════════════════════

def stage6_md(cfg: dict, paths: dict) -> dict:
    """
    Load AMBER prmtop into OpenMM, run:
      - Energy minimisation
      - NVT equilibration
      - NPT production MD

    UnoMD's force field parameterization step is intentionally bypassed;
    we inject the AMBER system directly into OpenMM.
    """
    log.info("=" * 60)
    log.info("STAGE 6 — OPENMM MD")
    log.info("=" * 60)

    omm = _import_openmm()
    AmberPrmtopFile  = omm["AmberPrmtopFile"]
    AmberInpcrdFile  = omm["AmberInpcrdFile"]
    DCDReporter      = omm["DCDReporter"]
    StateDataReporter = omm["StateDataReporter"]
    Simulation       = omm["Simulation"]
    LangevinMiddleIntegrator = omm["LangevinMiddleIntegrator"]
    MonteCarloBarostat = omm["MonteCarloBarostat"]
    Platform         = omm["Platform"]
    PME              = omm["PME"]
    HBonds           = omm["HBonds"]
    unit             = omm["unit"]

    wdir = Path(cfg["workdir"]) / "04_md"
    wdir.mkdir(parents=True, exist_ok=True)

    prmtop_file = paths["prmtop"]
    inpcrd_file = paths["inpcrd"]

    T    = cfg["temperature_K"] * unit.kelvin
    fric = cfg["friction"] / unit.picosecond
    dt   = cfg["timestep_fs"] * 0.001 * unit.picoseconds

    equil_steps = int(cfg["equil_ns"] * 1e3 / cfg["timestep_fs"] * 1000)
    prod_steps  = int(cfg["prod_ns"]  * 1e3 / cfg["timestep_fs"] * 1000)
    save_int    = cfg["save_interval"]

    log.info(f"  System:         {prmtop_file}")
    log.info(f"  Platform:       {cfg['platform']}")
    log.info(f"  Temperature:    {cfg['temperature_K']} K")
    log.info(f"  Timestep:       {cfg['timestep_fs']} fs")
    log.info(f"  Equilibration:  {cfg['equil_ns']} ns  ({equil_steps:,} steps)")
    log.info(f"  Production:     {cfg['prod_ns']} ns  ({prod_steps:,} steps)")

    # ── load topology ─────────────────────────────────────────────────────────
    log.info("  Loading AMBER topology into OpenMM...")
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(inpcrd_file)

    # ── build system ──────────────────────────────────────────────────────────
    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * unit.nanometers,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )

    # ── platform ─────────────────────────────────────────────────────────────
    try:
        platform = Platform.getPlatformByName(cfg["platform"])
        platform_props = {"Precision": cfg["precision"]} if cfg["platform"] == "CUDA" else {}
    except Exception:
        log.warning(f"  Platform {cfg['platform']} unavailable, falling back to CPU.")
        platform = Platform.getPlatformByName("CPU")
        platform_props = {}

    integrator = LangevinMiddleIntegrator(T, fric, dt)

    simulation = Simulation(prmtop.topology, system, integrator, platform, platform_props)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # ── energy minimisation ───────────────────────────────────────────────────
    log.info("  Energy minimisation...")
    tol = cfg["emin_tol"] * unit.kilojoule_per_mole / unit.nanometer
    simulation.minimizeEnergy(tolerance=tol, maxIterations=5000)
    emin_pdb = wdir / "emin.pdb"
    _save_pdb(simulation, prmtop, emin_pdb)
    log.info(f"  Emin structure → {emin_pdb}")

    # ── NVT equilibration ─────────────────────────────────────────────────────
    log.info(f"  NVT equilibration ({cfg['equil_ns']} ns)...")
    simulation.context.setVelocitiesToTemperature(T)

    equil_dcd = wdir / "equil.dcd"
    equil_log = wdir / "equil.log"
    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(str(equil_dcd), save_int))
    simulation.reporters.append(StateDataReporter(
        str(equil_log), save_int,
        step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, temperature=True, progress=True,
        totalSteps=equil_steps, speed=True,
    ))
    simulation.step(equil_steps)
    log.info(f"  Equilibration trajectory → {equil_dcd}")

    # ── NPT production ────────────────────────────────────────────────────────
    if cfg["npt_prod"]:
        log.info("  Adding Monte Carlo barostat for NPT production...")
        system.addForce(MonteCarloBarostat(1.0 * unit.bar, T, 25))
        simulation.context.reinitialize(preserveState=True)

    log.info(f"  Production MD ({cfg['prod_ns']} ns)...")
    prod_dcd = wdir / "production.dcd"
    prod_log = wdir / "production.log"
    chk_file = wdir / "production.chk"

    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(str(prod_dcd), save_int))
    simulation.reporters.append(StateDataReporter(
        str(prod_log), save_int,
        step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, temperature=True, density=True,
        progress=True, totalSteps=prod_steps, speed=True,
    ))
    # checkpoint every 1M steps (~2 ns)
    from openmm.app import CheckpointReporter
    simulation.reporters.append(CheckpointReporter(str(chk_file), 1000000))

    simulation.step(prod_steps)

    # save final state
    state_xml = wdir / "final_state.xml"
    simulation.saveState(str(state_xml))

    log.info(f"  Production trajectory → {prod_dcd}")
    log.info(f"  Final state           → {state_xml}")

    return {
        **paths,
        "md_wdir":   str(wdir),
        "emin_pdb":  str(emin_pdb),
        "equil_dcd": str(equil_dcd),
        "prod_dcd":  str(prod_dcd),
        "prod_log":  str(prod_log),
    }


def _save_pdb(simulation, prmtop, outpath: Path):
    from openmm.app import PDBFile
    state = simulation.context.getState(getPositions=True)
    with open(outpath, "w") as f:
        PDBFile.writeFile(prmtop.topology, state.getPositions(), f)


# ═════════════════════════════════════════════════════════════════════════════
#  STAGE 7 — ANALYSIS
# ═════════════════════════════════════════════════════════════════════════════

def stage7_analysis(cfg: dict, paths: dict) -> None:
    """
    PaDa-I-specific analysis:
      - Backbone RMSD vs. 5OXU reference
      - Per-residue RMSF
      - Fe–S bond length timeseries
      - Phe191 chi1 dihedral timeseries
      - Channel entrance width (Phe76–Phe191 CA distance)
    """
    log.info("=" * 60)
    log.info("STAGE 7 — ANALYSIS")
    log.info("=" * 60)

    mda, rms_mod, dih_mod, dist_fn = _import_mda()
    import numpy as np

    wdir = Path(cfg["workdir"]) / "05_analysis"
    wdir.mkdir(parents=True, exist_ok=True)

    prmtop  = paths["prmtop"]
    prod_dcd = paths["prod_dcd"]
    dry_pdb  = paths["dry_pdb"]

    log.info(f"  Loading trajectory: {prod_dcd}")
    u = mda.Universe(prmtop, prod_dcd)

    ca_sel     = u.select_atoms("protein and name CA")
    fe_sel     = u.select_atoms("resname HEM and name FE")
    sg_sel     = u.select_atoms(f"resid {cfg['cys_axial']} and name SG")
    phe191_sel = u.select_atoms("resid 191 and protein")
    phe76_ca   = u.select_atoms("resid 76  and name CA")
    phe191_ca  = u.select_atoms("resid 191 and name CA")

    # ── backbone RMSD ────────────────────────────────────────────────────────
    log.info("  Computing backbone RMSD...")
    ref = mda.Universe(dry_pdb)
    rmsd_analysis = rms_mod.RMSD(ca_sel, ref.select_atoms("protein and name CA"), select="backbone")
    rmsd_analysis.run()
    rmsd_out = wdir / "rmsd_backbone.csv"
    np.savetxt(rmsd_out, rmsd_analysis.results.rmsd[:, [1, 2]],
               delimiter=",", header="time_ps,rmsd_angstrom", comments="")
    log.info(f"  RMSD → {rmsd_out}")

    # ── per-residue RMSF ─────────────────────────────────────────────────────
    log.info("  Computing per-residue RMSF...")
    rmsf_analysis = rms_mod.RMSF(ca_sel)
    rmsf_analysis.run()
    rmsf_out = wdir / "rmsf_per_residue.csv"
    resids = ca_sel.resids
    with open(rmsf_out, "w") as f:
        f.write("resid,rmsf_angstrom\n")
        for res, rmsf_val in zip(resids, rmsf_analysis.results.rmsf):
            f.write(f"{res},{rmsf_val:.4f}\n")
    log.info(f"  RMSF → {rmsf_out}")

    # ── Fe–S bond length timeseries ──────────────────────────────────────────
    log.info("  Computing Fe–S bond length...")
    fe_s_distances = []
    times = []
    warn_thresh = cfg["fe_s_warn_ang"]
    n_warn = 0

    for ts in u.trajectory:
        if len(fe_sel) == 0 or len(sg_sel) == 0:
            log.warning("  Fe or SG atom not found — skipping Fe–S analysis.")
            break
        d = float(np.linalg.norm(fe_sel.positions[0] - sg_sel.positions[0]))
        fe_s_distances.append(d)
        times.append(ts.time)
        if d > warn_thresh:
            n_warn += 1

    if fe_s_distances:
        fe_s_out = wdir / "fe_s_bond.csv"
        np.savetxt(fe_s_out, np.column_stack([times, fe_s_distances]),
                   delimiter=",", header="time_ps,fe_s_angstrom", comments="")
        avg_d = np.mean(fe_s_distances)
        log.info(f"  Fe–S mean: {avg_d:.2f} Å  (frames > {warn_thresh} Å: {n_warn})")
        if n_warn > len(fe_s_distances) * 0.05:
            log.warning(
                f"  Fe–S bond drifted > {warn_thresh} Å in {n_warn} frames "
                f"({100*n_warn/len(fe_s_distances):.1f}%). "
                f"Check heme parameters."
            )
        log.info(f"  Fe–S timeseries → {fe_s_out}")

    # ── Phe191 chi1 dihedral (hinge residue) ────────────────────────────────
    log.info("  Computing Phe191 chi1 dihedral (N-CA-CB-CG)...")
    phe191_dih_sel = u.select_atoms("resid 191 and name N CA CB CG")
    if len(phe191_dih_sel) == 4:
        chi1_vals, chi1_times = [], []
        for ts in u.trajectory:
            # manual dihedral calculation
            p = phe191_dih_sel.positions
            b1 = p[1] - p[0]
            b2 = p[2] - p[1]
            b3 = p[3] - p[2]
            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)
            m1 = np.cross(n1, b2 / np.linalg.norm(b2))
            x = np.dot(n1, n2)
            y = np.dot(m1, n2)
            chi1_vals.append(np.degrees(np.arctan2(y, x)))
            chi1_times.append(ts.time)
        chi1_out = wdir / "phe191_chi1.csv"
        np.savetxt(chi1_out, np.column_stack([chi1_times, chi1_vals]),
                   delimiter=",", header="time_ps,chi1_degrees", comments="")
        log.info(f"  Phe191 chi1 → {chi1_out}")
    else:
        log.warning(f"  Phe191 dihedral atoms not found (got {len(phe191_dih_sel)}/4).")

    # ── channel entrance width (Phe76–Phe191 CA distance) ───────────────────
    log.info("  Computing channel entrance width (Phe76–Phe191 CA)...")
    if len(phe76_ca) == 1 and len(phe191_ca) == 1:
        cw_times, cw_dists = [], []
        for ts in u.trajectory:
            d = float(np.linalg.norm(phe76_ca.positions[0] - phe191_ca.positions[0]))
            cw_times.append(ts.time)
            cw_dists.append(d)
        cw_out = wdir / "channel_width_phe76_phe191.csv"
        np.savetxt(cw_out, np.column_stack([cw_times, cw_dists]),
                   delimiter=",", header="time_ps,width_angstrom", comments="")
        log.info(f"  Channel width mean: {np.mean(cw_dists):.2f} Å  "
                 f"(crystal ref: ~7.8 Å)")
        log.info(f"  Channel width → {cw_out}")
    else:
        log.warning("  Phe76 or Phe191 CA not found — skipping channel width.")

    log.info(f"  All analysis outputs in: {wdir}")


# ═════════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ═════════════════════════════════════════════════════════════════════════════

STAGES = {
    1: ("PDB preparation",         stage1_prepare_pdb),
    2: ("MCPB.py step 1",          stage2_mcpb_step1),
    3: ("Gaussian QM",             stage3_gaussian),
    4: ("MCPB.py steps 2–4",       stage4_mcpb_steps234),
    5: ("tleap system assembly",   stage5_tleap),
    6: ("OpenMM MD",               stage6_md),
    7: ("Analysis",                stage7_analysis),
}


def main():
    parser = argparse.ArgumentParser(
        description="Full PaDa-I UPO MD pipeline with MCPB.py heme parameterization.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Examples:
              # Run full pipeline (all stages):
              python pada1_md_pipeline.py

              # Run from stage 5 onwards (prmtop already built):
              python pada1_md_pipeline.py --stage 5

              # Use custom config and add a docked substrate:
              python pada1_md_pipeline.py --config my_config.yaml --ligand naphthalene.sdf

              # Analysis only:
              python pada1_md_pipeline.py --stage 7
        """)
    )
    parser.add_argument("--stage",  type=int, default=1,
                        help="Start from this stage (1–7). Default: 1")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to YAML config file (overrides defaults)")
    parser.add_argument("--ligand", type=str, default=None,
                        help="Path to substrate .sdf file (optional)")
    parser.add_argument("--workdir", type=str, default=None,
                        help="Working directory (overrides config)")
    parser.add_argument("--list-stages", action="store_true",
                        help="Print stage list and exit")
    args = parser.parse_args()

    if args.list_stages:
        print("\nAvailable stages:")
        for n, (name, _) in STAGES.items():
            print(f"  {n}: {name}")
        sys.exit(0)

    cfg = load_config(args.config)
    if args.ligand:
        cfg["ligand_sdf"] = args.ligand
    if args.workdir:
        cfg["workdir"] = args.workdir

    Path(cfg["workdir"]).mkdir(parents=True, exist_ok=True)

    log.info("╔══════════════════════════════════════════════════════════╗")
    log.info("║   PaDa-I UPO — Full MD Pipeline   BioBuilders 2026     ║")
    log.info("╚══════════════════════════════════════════════════════════╝")
    log.info(f"  Working directory: {cfg['workdir']}")
    log.info(f"  Starting stage:    {args.stage}")
    log.info(f"  Ligand:            {cfg.get('ligand_sdf') or 'none'}")

    # ── inter-stage state (paths dict passed through) ─────────────────────────
    paths = {}

    # if starting mid-pipeline, try to reconstruct paths from existing files
    if args.stage > 1:
        paths = _recover_paths(cfg, args.stage)

    for stage_num in range(args.stage, max(STAGES.keys()) + 1):
        name, fn = STAGES[stage_num]
        log.info(f"\n{'─'*60}")
        log.info(f"Running stage {stage_num}: {name}")
        log.info(f"{'─'*60}")

        if stage_num == 7:
            fn(cfg, paths)   # analysis returns None
        else:
            paths = fn(cfg, paths)

    log.info("\n✓ Pipeline complete.")
    log.info(f"  prmtop:       {paths.get('prmtop', 'n/a')}")
    log.info(f"  Production:   {paths.get('prod_dcd', 'n/a')}")
    log.info(f"  Analysis:     {cfg['workdir']}/05_analysis/")


def _recover_paths(cfg: dict, start_stage: int) -> dict:
    """
    Reconstruct the paths dict when resuming from a mid-pipeline stage.
    Looks for expected files in the workdir subdirectories.
    """
    pdb_id = cfg["pdb_id"].lower()
    wd     = Path(cfg["workdir"])
    paths  = {"pdb_id": pdb_id}

    candidates = {
        "prot_pdb":   wd / "01_prep" / f"{cfg['pdb_id']}_protein.pdb",
        "heme_pdb":   wd / "01_prep" / f"{cfg['pdb_id']}_heme.pdb",
        "mg_pdb":     wd / "01_prep" / f"{cfg['pdb_id']}_mg.pdb",
        "wdir":       wd / "02_mcpb",
        "mcpb_in":    wd / "02_mcpb" / f"{pdb_id}_heme.in",
        "small_opt":  wd / "02_mcpb" / f"{pdb_id}_heme_small_opt.com",
        "small_fc":   wd / "02_mcpb" / f"{pdb_id}_heme_small_fc.com",
        "large_mk":   wd / "02_mcpb" / f"{pdb_id}_heme_large_mk.com",
        "frcmod":     wd / "02_mcpb" / f"{pdb_id}_heme_small.frcmod",
        "resp_mol2":  wd / "02_mcpb" / f"{pdb_id}_heme_large_RESP.mol2",
        "tleap_in":   wd / "02_mcpb" / f"{pdb_id}_heme_tleap.in",
        "prmtop":     wd / "03_tleap" / f"{pdb_id}_solvated.prmtop",
        "inpcrd":     wd / "03_tleap" / f"{pdb_id}_solvated.inpcrd",
        "dry_pdb":    wd / "03_tleap" / f"{pdb_id}_dry.pdb",
        "prod_dcd":   wd / "04_md"    / "production.dcd",
        "prod_log":   wd / "04_md"    / "production.log",
        "emin_pdb":   wd / "04_md"    / "emin.pdb",
        "equil_dcd":  wd / "04_md"    / "equil.dcd",
        "md_wdir":    wd / "04_md",
    }

    for key, path in candidates.items():
        p = Path(path)
        if p.exists():
            paths[key] = str(p)
        else:
            log.debug(f"  Recovering paths: {key} not found at {path}")

    # convert wdir Path to str
    if "wdir" in paths:
        paths["wdir"] = str(paths["wdir"])

    required_for = {
        4: ["wdir", "mcpb_in", "small_opt", "small_fc", "large_mk"],
        5: ["prot_pdb", "frcmod", "resp_mol2"],
        6: ["prmtop", "inpcrd"],
        7: ["prmtop", "prod_dcd", "dry_pdb"],
    }
    if start_stage in required_for:
        missing = [k for k in required_for[start_stage] if k not in paths]
        if missing:
            log.warning(
                f"Resuming at stage {start_stage} but missing: {missing}\n"
                f"Run earlier stages first, or check workdir: {cfg['workdir']}"
            )

    return paths


if __name__ == "__main__":
    main()
