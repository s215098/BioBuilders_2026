# PaDa-I UPO — Full MD Pipeline

Automated molecular dynamics pipeline for **PaDa-I unspecific peroxygenase** (PDB: [5OXU](https://www.rcsb.org/structure/5OXU)), developed by DTU BioBuilders 2026 for the iGEM epoxy resin degradation project.

The pipeline covers everything from raw PDB download to production MD and analysis, including quantum-mechanical parameterization of the heme–Fe active site.

---

## Pipeline stages

| # | Stage | Tool |
|---|-------|------|
| 1 | PDB preparation — download, clean, protonate, split | `pdb4amber`, `pdbfixer` |
| 2 | MCPB.py step 1 — build QM model, generate GAMESS inputs | `MCPB.py` (AmberTools) |
| 3 | QM calculations — geometry opt + force constants + ESP charges | GAMESS-US |
| 4 | MCPB.py steps 2–4 — extract FF params, generate tleap input | `MCPB.py` |
| 5 | tleap — assemble solvated system → prmtop/inpcrd | `tleap` (AmberTools) |
| 6 | OpenMM — energy minimisation → NVT equilibration → production MD | OpenMM 8.x |
| 7 | Analysis — RMSD, RMSF, Fe–S bond, Phe191 dihedral, channel width | `MDAnalysis` |

---

## Files

```
MD/
├── pada1_md_pipeline.py      # main pipeline script
├── pada1_config.yaml         # production config (HPC, 100 ns)
├── pada1_test_config.yaml    # fast test config (1000 steps, CPU)
└── tests/
    └── test_pipeline_gamess.py   # unit tests for GAMESS/MCPB integration
```

---

## Requirements

### Conda environment

Create the `pada1_md` environment (requires AmberTools 24+ and OpenMM 8+):

```bash
conda create -n pada1_md -c conda-forge -c ambermd \
    python=3.11 ambertools=24 openmm=8 mdanalysis pdbfixer pyyaml
conda activate pada1_md
```

### GAMESS-US

Download and compile [GAMESS-US](https://www.msg.chem.iastate.edu/gamess/download.html) (free academic licence).  
The pipeline expects the `rungms` launcher — set `gamess_exe` in your config to its full path.

> GAMESS is **not included** in this repo due to licence restrictions.

---

## Quick start

### 1. Test run (CPU, ~minutes for stages 1–2)

```bash
conda activate pada1_md
python pada1_md_pipeline.py --config pada1_test_config.yaml --end-stage 2
```

### 2. Validate QM setup (~2–5 min single-point energy)

```bash
python pada1_md_pipeline.py --config pada1_test_config.yaml --validate-qm
```

This runs a GAMESS single-point on the small model to confirm the QM inputs are correct before submitting the full geometry optimisation.

### 3. Full pipeline

```bash
python pada1_md_pipeline.py --config pada1_config.yaml
```

### Resume from a specific stage

```bash
python pada1_md_pipeline.py --config pada1_config.yaml --stage 4
```

---

## Configuration

Key fields in `pada1_config.yaml` / `pada1_test_config.yaml`:

| Field | Description | Default |
|-------|-------------|---------|
| `pdb_id` | PDB entry to download | `5OXU` |
| `workdir` | Output directory | `pada1_md` |
| `qm_engine` | `gamess` or `gaussian` | `gamess` |
| `qm_method` | `pbe0_def2svp` (recommended) or `b3lyp_631g` | `pbe0_def2svp` |
| `gamess_exe` | Path to `rungms` | `/root/gamess/rungms` |
| `gamess_ncpus` | DDI process count | `1` (test) / `16` (HPC) |
| `platform` | `CUDA`, `CPU`, or `OpenCL` | `CUDA` |
| `prod_ns` | Production MD length (ns) | `100.0` |

### HPC (16-core) settings

See `pada1_config.yaml` — already configured for 16 DDI processes with 2000 MW local memory and 16000 MW distributed DDI pool (≈256 GB total).  
Update `gamess_exe` to the actual `rungms` path on your cluster.

---

## QM method: PBE0/SPK-DZP + SBK ECP on Fe

The three GAMESS jobs use:
- **PBE0** hybrid DFT (25 % HF exchange) — better than B3LYP for Fe spin-state energetics
- **SPK-DZP** Sapporo double-zeta polarization basis on all atoms (spherical harmonics, `ISPHER=1`)
- **SBK ECP** on Fe — replaces the Ar core (18 electrons), provides frozen-core treatment
- **ROHF** open-shell formalism for ferric high-spin Fe³⁺ (S = 5/2, multiplicity = 6)

---

## Tests

```bash
conda activate pada1_md
pip install pytest
pytest tests/test_pipeline_gamess.py -v
```

15 unit tests covering default config, MCPB.py input generation, GAMESS patching, and convergence detection.

---

## Output structure

```
pada1_md/
├── 01_prep/       — cleaned PDB, protein + heme + Mg split files
├── 02_mcpb/       — MCPB.py inputs, GAMESS .inp files, FF parameters
├── 03_tleap/      — solvated prmtop/inpcrd, dry PDB
├── 04_md/         — trajectory (.dcd), energy log, checkpoint PDBs
└── 05_analysis/   — RMSD/RMSF plots, Fe–S distance, dihedral timeseries
```
