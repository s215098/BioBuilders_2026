# Channel discovery pipeline

Two-stage structural screening for heme peroxygenases (UPOs). Given a set of
receptor structures, it (1) finds the heme-proximal active-site pocket and docks
a substrate into it, then (2) maps the substrate-access tunnels leading to that
pocket and compares their shape across structures.

The two stages live in their own folders and run **in order** — the tunnel stage
consumes the pocket stage's output:

```
channel_discovery/
├── fpocket_pipeline/   # STAGE A — pocket detection + docking
└── caver_pipeline/     # STAGE B — tunnel detection + clustering (needs Stage A output)
```

Each folder has its own detailed README:
[`fpocket_pipeline/README.md`](fpocket_pipeline/README.md) and
[`caver_pipeline/README.md`](caver_pipeline/README.md). This document is the
top-level map: what to run, in what order, and why.

---

## Why the order matters

`fpocket_pipeline` writes `outputs/pocket_selection.tsv` — a table listing, for
every structure, whether its active-site pocket passed quality checks (`PASS`/`FAIL`)
and which pocket number was selected. `caver_pipeline` reads that file via its
`--fpocket-dir` flag and:

- **skips** any structure that `FAIL`ed the pocket check (no tunnel analysis), and
- **filters** each structure's tunnels by proximity to the selected pocket.

So always run Stage A first, then point Stage B at its `outputs/` folder.

---

## Requirements

| Tool | Used by | Notes |
|------|---------|-------|
| Python ≥ 3.8 with `numpy`, `scipy`, `matplotlib` | both stages | pocket filtering, tunnel filtering, all plots |
| `rdkit`, `vina`, `meeko` | Stage A docking only | provided by the **`biobuilders`** conda env |
| `fpocket` | Stage A | `conda install -c bioconda fpocket` (Linux + macOS). A pre-built **Linux** binary is vendored at `fpocket_pipeline/fpocket/bin/fpocket` as a fallback |
| Java ≥ 1.8 | Stage B | runs `caver_pipeline/caver/caver.jar` |

The docking steps (Stage A, steps 3–4) require the conda environment:

```bash
conda activate biobuilders
```

> **fpocket is platform-specific.** The vendored binary is a Linux (ELF)
> executable — it will **not** run on macOS. The recommended cross-platform
> route is `conda install -c bioconda fpocket` (add it to the `biobuilders`
> env). The pipeline automatically uses an `fpocket` found on `PATH` and only
> falls back to the vendored Linux binary if none is installed, so Linux users
> can rely on the fallback while macOS users install via conda (or rebuild with
> `cd fpocket_pipeline/fpocket && make`).

Pocket detection/filtering (Stage A steps 1–2) and the entire tunnel stage
(Stage B) need only `numpy`/`scipy`/`matplotlib` + `fpocket`/`java`, so they can
be run without the conda environment.

---

## Quick start

Run everything from `channel_discovery/`. Put your receptor `.pdb`/`.cif` files in
`fpocket_pipeline/inputs/`. Three test structures (`5FUJ_heme`, `5OXU_heme`,
`9HE9_heme`) are already included.

```bash
# ── STAGE A — pockets + docking ──────────────────────────────────────────────
conda activate biobuilders          # needed for the docking steps
cd fpocket_pipeline
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --prepare
#   → writes outputs/pocket_selection.tsv and docking_out/*.pdbqt
cd ..

# ── STAGE B — tunnels + clustering ───────────────────────────────────────────
cd caver_pipeline
# copy your receptor PDBs here (or reuse the same structures)
cp ../fpocket_pipeline/inputs/*.pdb input_files/

python3 python_scripts/run_batch.py --fpocket-dir ../fpocket_pipeline/outputs
#   → runs CAVER on every PASS structure, writes results/tunnel_selection.csv

python3 python_scripts/cluster_tunnels.py --fpocket-dir ../fpocket_pipeline/outputs
#   → clusters tunnel shapes, writes results/clustering/cluster_assignments.csv
```

That is the full path from raw structures to clustered substrate-access channels.

---

## Stage A — pocket detection and docking (`fpocket_pipeline`)

Orchestrated by `pipeline.py`, which runs four steps in sequence:

1. **prepare** (`--prepare`) — convert each `.pdb`/`.cif` receptor to `.pdbqt`.
2. **fpocket** — detect pockets in every structure.
3. **filter** — for each structure, select the first pocket that is both
   heme-proximal and passes the metric thresholds.
4. **dock + analyze** — targeted AutoDock Vina docking into the selected pocket,
   then check each pose's orientation relative to the heme iron.

```bash
cd fpocket_pipeline
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --prepare
```

Useful variants:

```bash
# Only detect + filter pockets, no docking (no conda env needed)
python3 pipeline.py inputs/ --skip-docking

# Re-run without redoing fpocket (outputs/ already populated)
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --skip-fpocket
```

**Key output:** `outputs/pocket_selection.tsv` — the PASS/FAIL + selected-pocket
table that Stage B depends on. See
[`fpocket_pipeline/README.md`](fpocket_pipeline/README.md) for every stage flag,
the filter thresholds, and the results-table column reference.

---

## Stage B — tunnel detection and analysis (`caver_pipeline`)

Run after Stage A. Point every command at the Stage A outputs with
`--fpocket-dir ../fpocket_pipeline/outputs`.

**1. Detect and filter tunnels** — runs CAVER on each PASS structure, injects the
heme Fe as the tunnel starting point, and keeps only tunnels close to the
selected pocket:

```bash
cd caver_pipeline
python3 python_scripts/run_batch.py --fpocket-dir ../fpocket_pipeline/outputs
#   → results/<stem>/ and results/tunnel_selection.csv
```

**2. Cluster tunnel shapes across structures** — groups structures by the shape
profile of their best tunnel:

```bash
python3 python_scripts/cluster_tunnels.py --fpocket-dir ../fpocket_pipeline/outputs
#   → results/clustering/cluster_assignments.csv (+ dendrograms and overlays)
```

**3. (Optional) Downstream analyses** — these read `cluster_assignments.csv`, so
run them after step 2:

```bash
# Find a discriminating stretch of the channel to cluster on, then re-cluster:
python3 python_scripts/find_reference_region.py
python3 python_scripts/cluster_tunnels.py --region <LO> <HI> --out results/clustering_region

# Link predicted binding affinity to tunnel clusters (needs a Boltz-2 summary CSV):
python3 python_scripts/cluster_affinity_violin.py --summary-csv <boltz_summary.csv>

# Compare tunnel clusters against CUPP sequence clades (needs a CUPP clade file):
python3 python_scripts/cluster_clade_heatmap.py
```

See [`caver_pipeline/README.md`](caver_pipeline/README.md) for the CAVER
parameters, the tunnel filters, clustering options, and region-restricted
clustering.

---

## Configuration

- **Stage A filter thresholds** (druggability, volume, Fe-proximity, …) are
  `pipeline.py` flags — see the fpocket README.
- **Stage B CAVER parameters** (probe radius, clustering threshold, …) live in
  `caver_pipeline/config.txt`, the shared template applied to every structure.
  Its `starting_point_residue` line is overwritten automatically per structure.

---

## Test data

The three included structures — `5FUJ_heme`, `5OXU_heme`, `9HE9_heme` — are enough
to exercise both stages end to end. Replace them with your own receptors to run a
real screen.
