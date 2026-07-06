# Pocket-guided docking for heme peroxygenases (UPOs)

End-to-end screening pipeline: given a folder of receptor PDB or mmCIF files,
it runs fpocket, selects the heme-proximal pocket for each structure, docks a
substrate with AutoDock Vina, and checks whether each pose has the correct
orientation relative to the heme iron.

Currently configured for **NNBT** as the substrate and **8AV5** as the
reference structure (druggability score 0.812, active-site pocket sitting
directly above the heme Fe in the aromatic "Phe cage" typical of UPO active
sites).

---

## Directory layout

```
inputs/                  # receptor PDB/CIF files and their PDBQTs
  8AV5_heme.pdb
  8AV5_heme.pdbqt
ligands/                 # ligand SDF files
  NNBT.sdf
outputs/                 # fpocket results, one sub-folder per receptor
  8AV5_heme_out/
    8AV5_heme_info.txt   #   per-pocket scores and descriptors
    pockets/
      pocket1_vert.pqr   #   Pocket 1 alpha spheres
      pocket1_atm.pdb    #   Pocket 1 lining atoms
      ...
  pocket_selection.tsv   # full results table written after every pipeline run
docking_out/             # Vina docking poses (PDBQT, one file per run)
pipeline.py              # full end-to-end orchestrator
prepare_receptor.py      # PDB/CIF → PDBQT conversion (single file or batch)
run_fpocket.py           # batch fpocket runner (PDB and mmCIF)
filter_pockets.py        # pocket metric filtering; heme-proximal selection when called via pipeline.py
pocket_box.py            # compute Vina search box from a pocket
pocket_residues.py       # list residues lining a pocket
docking.py               # run AutoDock Vina
analyze_poses.py         # check docked poses for feasibility and orientation
fpocket/                 # fpocket tool source and compiled binaries
```

---

## mmCIF support

The pipeline accepts `.cif` files alongside `.pdb` files throughout:

- **`prepare_receptor.py`** — auto-detects mmCIF space-separated format (as
  produced by RFdiffusion/ProteinMPNN) and converts to PDBQT. Standard PDB
  column format is handled the same way.
- **`run_fpocket.py`** — for CIF inputs, fpocket is run on the PDBQT (which is
  PDB-column-compatible) via a temp file with a sanitised name. This works
  around fpocket's format detector, which incorrectly routes any filename
  containing `.cif` to its CIF parser. The output folder is named after the
  original CIF stem as normal.
- **`pipeline.py` `--prepare`** — converts both `.pdb` and `.cif` inputs.
- **Heme Fe detection for CIF inputs** — `find_heme_fe()` reads PDB columns, so
  for CIF inputs the pipeline falls back to fpocket's own `<stem>_out.pdb`,
  which fpocket always writes in PDB format regardless of input type.

---

## Requirements

All scripts require the **`biobuilders`** conda environment:

```bash
conda activate biobuilders
```

`pocket_box.py`, `pocket_residues.py`, `run_fpocket.py`, `filter_pockets.py`,
and `prepare_receptor.py` need only `numpy`. `docking.py` and `analyze_poses.py`
additionally require `vina`, `rdkit`, and `meeko`.

---

## Quick start — full pipeline

```bash
source /home/touliopoulos/miniconda3/etc/profile.d/conda.sh
conda activate biobuilders
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --prepare
```

`--prepare` converts any receptor without a matching `.pdbqt` before running.
Drop it if all PDBQTs already exist.

---

## Pipeline stages

### Stage 0 — `prepare_receptor.py` (opt-in via `--prepare`)

Converts PDB or mmCIF → PDBQT for any receptor that lacks a matching `.pdbqt`
in `input_dir`. Handles ATOM and HETATM records (including heme and any HETATM
named FE), removes crystallographic waters by default, and strips alternate
conformations. Auto-detects format (PDB column layout vs. mmCIF space-separated).

```bash
# Single file (PDB or CIF)
python3 prepare_receptor.py -i inputs/MyUPO.pdb
python3 prepare_receptor.py -i inputs/design.cif

# Batch (all .pdb and .cif files in a directory)
python3 prepare_receptor.py -i inputs/

# Keep crystallographic waters
python3 prepare_receptor.py -i inputs/ --keep-water
```

| Option | Meaning | Default |
|--------|---------|---------|
| `-i` | PDB/CIF file or directory | required |
| `-o` | Output PDBQT or directory | same location as input |
| `--keep-water` | Retain HOH/WAT molecules | off |
| `--overwrite` | Re-convert existing PDBQTs | off |

---

### Stage 1 — `run_fpocket.py`

Runs fpocket on every `.pdb` and `.cif` file in a directory. Each `<stem>_out/`
folder is saved under `outputs/`. For CIF-derived inputs, uses the pre-built
PDBQT via a temp file to avoid fpocket's filename-based format detection bug.

```bash
python3 run_fpocket.py inputs/
python3 run_fpocket.py inputs/ -o outputs/ --overwrite
```

| Option | Meaning | Default |
|--------|---------|---------|
| `input_dir` | Folder containing PDB/CIF files | required |
| `-o, --output_dir` | Where to save `_out/` folders | `outputs/` |
| `--fpocket` | Path to fpocket binary | `fpocket/bin/fpocket` |
| `--overwrite` | Re-run even if output exists | off |

---

### Stage 2 — `filter_pockets.py` (heme-proximal pocket selection)

For each structure, iterates through all fpocket pockets in ranked order and
selects the **first pocket that is both heme-proximal and passes all metric
thresholds**. This means the pipeline no longer blindly uses Pocket 1 — if the
active-site pocket is ranked 2nd or 3rd by fpocket score, it is still found and
used for docking.

#### Pocket selection logic

```
for pocket_n in 1, 2, 3, ... (fpocket rank order):
    if min_distance(alpha_spheres → heme_Fe) > --fe-proximity  → skip
    if total_sasa == 0                                          → skip (buried)
    if any metric threshold fails                               → skip
    else → use this pocket for docking box; stop
```

The selected pocket number is reported in the Stage 2 output and propagates to
Stage 3 for box computation.

#### Heme iron detection

`find_heme_fe()` searches in this order:
1. Fe atom in a standard heme residue (HEM / HEC / HEB)
2. Centroid of all atoms in a standard heme residue
3. Any Fe atom in any HETATM record (catches non-standard names such as `LIG`)

If no Fe is found, the structure is skipped with a warning.

#### Metric thresholds

| Flag | Metric | Default | Why |
|------|--------|---------|-----|
| `--fe-proximity` | Max distance (Å) from any alpha sphere to heme Fe | ≤ 10.0 | Identifies the active-site access channel; non-active-site pockets are typically > 15 Å from Fe |
| `--druggability` | Druggability score (0–1) | ≥ 0.2 | ML estimate of drug-bindability |
| `--volume` | Pocket volume (Å³) | ≥ 500 | Practical minimum for a UPO active site; well above NNBT molecular volume (~187 Å³) |
| `--total-sasa` | Total pocket SASA (Å²) | ≥ 80 | Guards against buried cavities with no solvent access |
| `--apolar-frac` | Apolar alpha-sphere fraction | ≥ 0.4 | UPO active sites are hydrophobic (Phe cage) |
| `--n-spheres` | Number of alpha spheres | ≥ 15 | Minimum pocket definition |

For the 8AV5 reference, the active-site Pocket 1 has alpha spheres at **3.6 Å**
from Fe; the next closest pocket is at 11.9 Å. A cutoff of 10 Å cleanly
separates active-site from non-active-site pockets with margin for structural
variation. Use `--fe-proximity 5` for a stricter selection if your structures
have well-converged active-site geometry.

The heme-proximal pocket selection logic runs only when invoked via
`pipeline.py`. The standalone CLI applies metric filters to **Pocket 1** of
each structure, which is useful for a quick metric audit of your fpocket
outputs but does not iterate through pockets or check Fe proximity.

```bash
# Standalone use — applies metric filters to Pocket 1 of all outputs
python3 filter_pockets.py
python3 filter_pockets.py outputs/ --druggability 0.4 --volume 500
python3 filter_pockets.py --out-tsv metrics.tsv --out-list passing.txt
```

---

### Stage 3 — `docking.py`

Targeted AutoDock Vina docking. The pipeline auto-computes the search box from
the selected pocket's `_vert.pqr` with heme-distance clipping. No manual box
entry is needed.

```bash
# Manual targeted run (use pocket_box.py to get center/size first)
python3 docking.py \
    -r inputs/8AV5_heme.pdbqt \
    -l ligands/NNBT.sdf \
    -o docking_out \
    --center -38.41 21.99 -15.28 \
    --box_size 15.7 13.4 13.0 \
    --n_poses 10
```

| Option | Meaning | Default |
|--------|---------|---------|
| `-r, --receptor` | Receptor PDBQT | required |
| `-l, --ligand` | Ligand SDF | required |
| `-o, --out_dir` | Output directory | required |
| `--center X Y Z` | Box center; omit for blind docking | none |
| `--box_size X Y Z` | Box dimensions (Å) | `15 15 15` |
| `--exhaustiveness` | Search effort | `64` |
| `--n_poses` | Poses to write | `20` |

#### Manual docking for a specific pocket

If you want to dock against a pocket other than the one the pipeline selected,
compute its box manually and run docking directly:

```bash
# Step 1 — get box parameters for pocket 2
python3 pocket_box.py outputs/MyUPO_out/pockets/pocket2_vert.pqr \
    --heme-pdb inputs/MyUPO.pdb --max-dist 12

# Step 2 — dock with those parameters
python3 docking.py -r inputs/MyUPO.pdbqt -l ligands/NNBT.sdf \
    -o docking_out --center X Y Z --box_size SX SY SZ

# Step 3 — check orientation
python3 analyze_poses.py docking_out/MyUPO_NNBT_targeted_*.pdbqt \
    --receptor inputs/MyUPO.pdb
```

---

### Stage 4 — `analyze_poses.py`

Checks each Vina pose for two criteria specific to UPO substrate binding.

```bash
python3 analyze_poses.py docking_out/*.pdbqt --receptor inputs/8AV5_heme.pdb
```

#### How the checks work

**Check 1 — Feasibility**

```
d_ring  ≤  --max-ring-dist   (default 8.0 Å)
```

`d_ring` is the distance from the aromatic ring centroid (PDBQT atom type `A`)
to the heme Fe. Rejects poses too far from the catalytic iron to be oxidised.

**Check 2 — Orientation**

```
d_ring  <  d_oxy + --flexibility   (default flexibility 2.0 Å)
```

`d_oxy` is the distance from the hydroxyl oxygen centroid (type `OA`) to Fe.
Enforces that the **ring faces the heme** and the **hydroxyls face away** — the
geometry required for UPO-mediated oxidation. `flexibility` is a tolerance
allowing slight tilting without hard failure.

**Score gate**

```
Vina score  ≤  --min-score   (default -5.0 kcal/mol)
```

Applied before geometric checks. A pose passes only if **all three** conditions
hold simultaneously.

| Option | Meaning | Default |
|--------|---------|---------|
| `--receptor` | PDB/PDBQT with heme (for Fe location) | required |
| `--min-score` | Max acceptable Vina score (kcal/mol) | `-5.0` |
| `--max-ring-dist` | Max ring-to-Fe distance (Å) | `8.0` |
| `--flexibility` | Orientation tolerance (Å) | `2.0` |
| `--top-n` | Only evaluate top N poses | all |

---

## `pipeline.py` — orchestrator

Runs all four stages in sequence. Existing docking outputs are skipped
automatically so re-runs only process new structures. At the end of every run
(including `--skip-docking` runs) a results table is written to
`outputs/pocket_selection.tsv`.

```bash
# Full pipeline with receptor preparation
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --prepare

# Skip fpocket (outputs/ already populated), tighten filter
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --skip-fpocket \
    --druggability 0.4 --volume 500

# Filtering only — no docking
python3 pipeline.py inputs/ --skip-fpocket --skip-docking

# Stricter heme-proximity cutoff and orientation check
python3 pipeline.py inputs/ -l ligands/NNBT.sdf --skip-fpocket \
    --fe-proximity 5 --min-score -6.0 --max-ring-dist 7.0
```

| Option group | Flags |
|---|---|
| Stages | `--prepare`, `--skip-fpocket`, `--skip-filter`, `--skip-docking` |
| Directories | `--outputs-dir`, `--docking-out` |
| Box | `--heme-max-dist` (default 12), `--buffer` (default 6), `--box-size X Y Z` |
| Vina | `--exhaustiveness` (default 32), `--n-poses` (default 10) |
| Filter | `--fe-proximity`, `--druggability`, `--volume`, `--total-sasa`, `--apolar-frac`, `--n-spheres` |
| Orientation | `--min-score`, `--max-ring-dist`, `--flexibility` |

---

## Results table — `outputs/pocket_selection.tsv`

Written at the end of every pipeline run, covering **all** structures regardless
of where they dropped out. Open in Excel, LibreOffice Calc, or:

```bash
column -t -s $'\t' outputs/pocket_selection.tsv | less -S
```

| Column | Meaning |
|--------|---------|
| `stem` | Structure name (filename without extension) |
| `filter_status` | `PASS` or `FAIL` from Stage 2 |
| `pocket` | Pocket number selected as heme-proximal (empty if FAIL) |
| `druggability` | Druggability score of the selected/candidate pocket |
| `volume` | Pocket volume (Å³) |
| `total_sasa` | Total solvent-accessible surface area (Å²) |
| `apolar_frac` | Apolar alpha-sphere fraction |
| `n_spheres` | Number of alpha spheres |
| `filter_reason` | Which threshold failed and on which pocket (empty if PASS) |
| `docking_status` | `already_docked`, `docked`, `no_pdbqt`, `box_failed`, `docking_failed`, `not_reached`, `skipped` |
| `n_poses` | Total Vina poses evaluated |
| `n_passing_poses` | Poses passing score + orientation criteria |
| `best_score` | Vina score of the best passing pose (kcal/mol) |
| `best_d_ring` | Ring centroid → Fe distance of the best passing pose (Å) |
| `best_d_oxy` | Hydroxyl centroid → Fe distance of the best passing pose (Å) |
| `outcome` | Final verdict — see table below |

**`outcome` values:**

| Value | Meaning |
|-------|---------|
| `pass` | ≥ 1 pose passed all score + orientation criteria |
| `no_passing_poses` | Docked successfully but no pose passed |
| `filtered_out` | Failed Stage 2 pocket filter |
| `no_pdbqt` | Passed filter but no `.pdbqt` found for docking |
| `box_failed` | Box computation failed (pocket PQR missing or heme issue) |
| `docking_failed` | Vina returned an error |
| `skipped` | Stage 3 not reached (`--skip-docking`) |

When `filter_reason` is populated for a FAIL row, it includes the candidate
pocket number, its distance to Fe, and the specific metric(s) that failed —
e.g. `pocket 3 (min_dist=3.4Å): volume 588.2 < 600.0`. This lets you quickly
judge whether to relax a threshold for a borderline structure.

---

## Utility scripts

### `pocket_box.py` — compute the docking box manually

```bash
python3 pocket_box.py outputs/8AV5_heme_out/pockets/pocket1_vert.pqr \
    --heme-pdb inputs/8AV5_heme.pdb --max-dist 12
```

### `pocket_residues.py` — list residues lining a pocket

```bash
python3 pocket_residues.py outputs/8AV5_heme_out/pockets/pocket1_atm.pdb

# Only the aromatic cage
python3 pocket_residues.py outputs/8AV5_heme_out/pockets/pocket1_atm.pdb \
    --resname PHE --resname TYR
```
