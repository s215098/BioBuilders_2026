# CAVER Pipeline

Tunnel analysis of peroxidase structures using [CAVER 3.0](http://www.caver.cz/).

---

## Directory structure

```
caver/
├── caver/               # CAVER tool (jar + binaries) — do not modify
├── config.txt           # Shared parameter template (probe radius, clustering, etc.)
├── input_files/         # Put your .pdb files here
├── results/             # Auto-created; one subfolder per PDB
│   ├── <pdb_stem>/
│   │   ├── config.txt        # Auto-generated config pointing to the Fe starting residue
│   │   ├── pdb/<stem>.pdb    # Modified PDB with synthetic Fe ATOM record injected
│   │   ├── out/              # Raw CAVER output (CSVs, PDB data, visualisation files)
│   │   └── plots/            # Generated PNG plots (filtered tunnels only)
│   ├── clustering/      # Cross-structure clustering plots
│   └── tunnel_selection.csv  # Per-tunnel filter summary across all structures
├── python_scripts/
│   ├── run_batch.py         # Batch runner: CAVER + filtering + plots for every PDB
│   ├── extract_tunnel_R.py  # Single-tunnel analysis and plotting
│   ├── cluster_tunnels.py   # Cross-structure clustering of tunnel shape profiles
│   ├── find_reference_region.py    # Locate the region where reference channels differ most
│   ├── cluster_affinity_violin.py  # Violin plots of predicted affinity across clusters
│   └── cluster_clade_heatmap.py    # Heatmap of tunnel clusters vs CUPP sequence clades
└── out/                 # Output of the original single-run (kept for reference)
```

---

## Requirements

- **Java** ≥ 1.8
- **Python** ≥ 3.8 with `matplotlib`, `numpy`, `scipy`

---

## Configuration

`config.txt` at the root is the **parameter template** used by every run.  
Edit it to change probe radius, clustering threshold, and other CAVER settings.

```
probe_radius          0.4
shell_radius          5
shell_depth           6
clustering_threshold  3.5
...
```

The `starting_point_residue` line is **overwritten automatically** for each PDB — you do
not need to set it manually. See [Starting point detection](#starting-point-detection) below.

---

## Running CAVER

### Single PDB (manual)

Run from the project root (`caver/`):

```bash
java -Xmx1200m -cp caver/lib -jar caver/caver.jar \
  -home caver/ \
  -pdb ./input_files/MyProtein.pdb \
  -conf ./config.txt \
  -out ./out
```

Output lands in `./out/`. Note that this uses the shared `config.txt` as-is, so
make sure `starting_point_residue` is set correctly for your structure before running.

---

### Batch — all PDBs in a folder

Place all `.pdb` files in `input_files/` (or any folder), then run from the project root:

```bash
python3 python_scripts/run_batch.py
```

For each PDB the script will:
1. Locate the Fe atom in the HEME/LIG group
2. Write a modified PDB (`results/<stem>/pdb/<stem>.pdb`) with the Fe injected as a synthetic `ATOM` record so CAVER can use its exact coordinates as the starting point
3. Generate `results/<stem>/config.txt` pointing `starting_point_residue` to that synthetic record
4. Run CAVER → output saved to `results/<stem>/out/`
5. Filter tunnels (see [Tunnel filtering](#tunnel-filtering) below)
6. Generate plots for kept tunnels → saved to `results/<stem>/plots/`
7. Write `results/tunnel_selection.csv` summarising every tunnel across all structures

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--pdb PATH` | — | Run on a single PDB file instead of a whole folder |
| `--pdb-dir PATH` | `input_files/` | Folder containing `.pdb` files |
| `--config PATH` | `config.txt` | Parameter template |
| `--results-dir PATH` | `results/` | Root folder for all outputs |
| `--fpocket-dir PATH` | — | fpocket outputs folder containing `pocket_selection.tsv` |
| `--pocket-radius FLOAT` | `12.0` | Max **mean** distance (Å) of a tunnel's path from the fpocket pocket center; tunnels above this are discarded |
| `--max-curvature FLOAT` | `1.5` | Discard tunnels whose curvature exceeds this value; set to a large number to disable |
| `--max-length FLOAT` | `25.0` | Discard tunnels longer than this value (Å); set to a large number to disable |
| `--probe-radius FLOAT` | template value | Override probe radius for CAVER (e.g. `0.4` to find tighter tunnels) |
| `--ylim FLOAT` | `5` | Y-axis half-range (±) for tunnel shape plots |
| `--skip-existing` | off | Skip both the CAVER run **and** plot regeneration for structures that already have output (`results/<stem>/out/summary.txt`); their rows are still added to the combined CSV. Use to add new PDBs without recomputing or re-plotting existing ones |
| `--no-caver` | off | Skip CAVER, regenerate plots and CSV only |
| `--no-plots` | off | Skip plot generation after CAVER runs |

#### Examples

```bash
# All PDBs, guided by fpocket pocket selection (curvature filter 1.5 applied by default)
python3 python_scripts/run_batch.py \
  --fpocket-dir ../fpocket_pipeline/outputs

# Relax curvature threshold to also include more curved paths
python3 python_scripts/run_batch.py \
  --fpocket-dir ../fpocket_pipeline/outputs \
  --max-curvature 2.0

# Disable curvature filter entirely
python3 python_scripts/run_batch.py \
  --fpocket-dir ../fpocket_pipeline/outputs \
  --max-curvature 99

# Regenerate plots/CSV from existing CAVER output (no re-run)
python3 python_scripts/run_batch.py \
  --no-caver \
  --fpocket-dir ../fpocket_pipeline/outputs

# Incremental run — after adding new PDBs (and/or updating fpocket): compute only the
# new PASS structures, reuse existing CAVER output, and rebuild the combined CSV
python3 python_scripts/run_batch.py \
  --fpocket-dir ../fpocket_pipeline/outputs \
  --skip-existing

# Use the original looser probe radius if 0.4 finds too many noise tunnels
python3 python_scripts/run_batch.py \
  --pdb input_files/MyProtein.pdb \
  --probe-radius 0.7
```

#### Key CAVER parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `probe_radius` | `0.4` Å | Minimum bottleneck radius a tunnel must have to be detected. Lower finds narrower tunnels but more noise. The default of 0.4 Å was chosen to detect tunnels in tightly packed AlphaFold models that produce no results at the commonly used 0.7 Å. |
| `clustering_threshold` | `3.5` Å | How different two tunnels must be to be reported as separate clusters. Raise to merge similar paths and reduce the number reported. |
| `shell_radius` | `5.0` Å | Search radius around each atom for Voronoi vertices. Rarely needs changing. |
| `shell_depth` | `6.0` Å | How far from the starting sphere CAVER extends its search. Rarely needs changing. |

---

## Starting point detection

CAVER requires a starting residue whose center of mass defines the center of the initial
search sphere. For substrate access tunnels in peroxidases, this should be the heme iron.

The script locates the Fe atom in any `HETATM` record belonging to a recognized heme group
(`HEM`, `HEC`, `LIG`, etc.) and writes a modified copy of the PDB with a synthetic `ATOM`
record placed at the exact Fe coordinates. `starting_point_residue` in the generated config
always points to this record, regardless of whether the structure is a crystal structure
(HEM in chain A) or an AlphaFold model (LIG in chain B). This ensures CAVER's starting
sphere is centered precisely on the iron for every structure.

---

## Tunnel filtering

CAVER often finds tens of tunnels per structure, most of which are irrelevant detours through
the protein far from the active site. Two sequential filters reduce this to the biologically
meaningful substrate access channels.

### Stage 1 — Structure-level (fpocket PASS/FAIL)

Requires `--fpocket-dir`. Structures that failed fpocket pocket quality checks (e.g. too
buried, too small) are **skipped entirely** — no CAVER run, no plots, no CSV rows. The
allow-list is read from `pocket_selection.tsv` inside the fpocket outputs folder.

### Stage 2 — Tunnel-level (proximity to the active-site pocket)

For each PASS structure, every tunnel found by CAVER is evaluated against the corresponding
fpocket pocket:

**Pocket-distance filter** (`--pocket-radius`, default 12 Å)  
The mean 3D distance of all points along the tunnel's path to the fpocket pocket centroid is
computed. Tunnels with a mean distance above the threshold are discarded. Using the mean
(rather than the minimum) prevents tunnels that start near the iron but veer off into
unrelated parts of the protein from passing — a tunnel that spends most of its length far
from the pocket will have a high mean even if its first point is close.

**Curvature filter** (`--max-curvature`, default `1.5`)  
CAVER's curvature metric is the ratio of tunnel length to the straight-line distance between
endpoints (1.0 = perfectly straight, 2.0 = twice as long as a direct path). Direct substrate
access tunnels typically have curvature 1.2–1.4; winding detour paths are often 1.8–3.0+.
The default of 1.5 keeps only the straighter, more direct channels. Moderately curved
(ellipsoidal) paths that stay within the pocket region still pass because their mean pocket
distance remains low. To disable this filter, pass a large value (e.g. `--max-curvature 99`).

**Length filter** (`--max-length`, default `25` Å)  
Discards tunnels longer than 25 Å. Very long tunnels are typically winding paths through
unrelated parts of the protein rather than direct substrate access channels. Note that
AlphaFold models tend to produce slightly longer tunnels than crystal structures (the iron
sits deeper in the predicted fold), so a threshold of 25 Å accommodates both.

A tunnel is **kept** only if it passes all three active filters. The `tunnel_selection.csv` records
the result for every tunnel with columns:

| Column | Description |
|--------|-------------|
| `stem` | Structure identifier |
| `fpocket_status` | `PASS` or empty (no fpocket data) |
| `pocket_num` | fpocket pocket number used for filtering |
| `tunnel` | CAVER tunnel ID |
| `filter_status` | `kept` or `discarded` |
| `discard_reason` | Which filter(s) triggered and the measured value |
| `throughput` | CAVER throughput score (higher = more open path) |
| `bottleneck_radius` | Narrowest radius along the tunnel (Å) |
| `length` | Tunnel length (Å) |
| `curvature` | Length / straight-line distance ratio |
| `pocket_dist` | Mean distance of tunnel path to pocket centroid (Å) |

Rows are ordered: kept tunnels first (highest throughput first), then discarded tunnels.

---

## Inspecting a single tunnel

`extract_tunnel_R.py` reads a `tunnel_profiles.csv` and produces plots for one tunnel:

```bash
python3 python_scripts/extract_tunnel_R.py --tunnel 1
```

| Flag | Default | Description |
|------|---------|-------------|
| `--tunnel INT` | `1` | Tunnel ID to analyse |
| `--csv PATH` | `out/analysis/tunnel_profiles.csv` | Path to tunnel_profiles.csv |
| `--out PATH` | `python_scripts/` | Output directory for plots |
| `--ylim FLOAT` | `3` | Y-axis half-range for the shape plot |

Two PNG files are produced per tunnel:

- `tunnel_<N>_R_distribution.png` — histogram of R values (frequency distribution)
- `tunnel_<N>_R_shape.png` — tunnel cross-section profile along distance

---

## Cross-structure tunnel clustering

`cluster_tunnels.py` compares the R-value shape profile of the best kept tunnel from each
structure and groups them by similarity using hierarchical clustering. Run it **after**
`run_batch.py` has produced CAVER output and `tunnel_selection.csv`.

**Tunnel selection for clustering:** the script reads `results/tunnel_selection.csv` and
picks the first `kept` tunnel per structure — i.e. the highest-throughput tunnel that passed
all active filters. Each structure may therefore contribute a different tunnel ID. If
`tunnel_selection.csv` is not present it falls back to tunnel 1. Use `--tunnel N` to force
a specific ID for all structures.

**Typical workflow:**

```bash
# Step 1 — run CAVER and filtering on all structures
python3 python_scripts/run_batch.py \
  --fpocket-dir ../fpocket_pipeline/outputs

# Step 2 — cluster the best kept tunnel across PASS structures
python3 python_scripts/cluster_tunnels.py \
  --fpocket-dir ../fpocket_pipeline/outputs

# Step 3 — inspect best_tunnel_dendrogram.png, choose a cluster count, rerun
python3 python_scripts/cluster_tunnels.py \
  --fpocket-dir ../fpocket_pipeline/outputs \
  --n-clusters 3
```

Each tunnel profile is resampled to a fixed-length normalized distance axis [0 = heme, 1 = exit]
so tunnels of different lengths can be directly compared. The entry radius, exit radius, and
bottleneck radius (minimum R) are appended as explicit features so they are never lost to
grid discretization. When `--fpocket-dir` is given, only fpocket PASS structures are included.

#### How profiles are normalized and resampled

CAVER measures each tunnel at an **uneven number of points at uneven spacings** — one structure's
tunnel might have 40 measurements, another's 90. Before they can be compared or clustered, every
profile is put on a common footing in two steps (see `build_feature_vector` in `cluster_tunnels.py`).

**Step 1 — normalize distance to 0–1.** CAVER reports a cumulative distance along the path (0 Å at
the heme, up to the full length at the exit). Dividing every distance by the total length rescales
the axis so all tunnels run from `0.0` (heme) to `1.0` (exit), regardless of their real length in Å.
Now "0.5" means "halfway along this particular tunnel" for every structure.

**Step 2 — resample onto the same fixed positions.** We pick a fixed set of evenly-spaced positions
(100 by default, `--n-points`) and ask every tunnel for its radius at those same positions. Where a
tunnel wasn't measured at exactly one of those spots, the value is **interpolated** — i.e. estimated
by connecting the two nearest real measurements with a straight line ("connect the dots") and reading
off the value in between.

A worked example with just 5 target positions (the code uses 100):

```
Two tunnels, measured at different points (already normalized to 0–1):

  Tunnel A:  position  0.0   0.33   0.66   1.0
             radius    1.0   2.5    3.0    2.0

  Tunnel B:  position  0.0   0.2   0.4   0.6   0.8   1.0
             radius    1.2   2.0   2.8   3.1   2.9   2.2

Read each tunnel's radius at the SAME 5 positions (0.0, 0.25, 0.5, 0.75, 1.0),
interpolating where there was no direct measurement:

             position  0.0   0.25   0.5   0.75   1.0
  Tunnel A:  radius    1.0   2.1    2.7   2.7    2.0
  Tunnel B:  radius    1.2   2.3    2.9   3.0    2.2
```

Both tunnels are now described by the same number of values at the same positions, so they can be
compared point-by-point and fed into clustering.

**Endpoints stay exact.** Because the target positions include `0.0` and `1.0` — which are real
measured points (the heme entrance and the exit) — no interpolation happens at the ends. Each
tunnel's true entrance and exit radii pass through unchanged. As a final safeguard the entry, exit,
and bottleneck (minimum) radii are also appended verbatim to the feature vector, so a narrow point
sitting between two grid positions can never be smoothed away before clustering.

The following files are saved to `results/clustering/`:

| File | Description |
|------|-------------|
| `best_tunnel_dendrogram.png` | Horizontal hierarchical tree — labels show stem and tunnel ID; sized to be square and readable for any number of structures |
| `best_tunnel_overlay.png` | All profiles overlaid and coloured by cluster |
| `best_tunnel_overlay_uncolored.png` | All profiles overlaid in a single colour, no clustering (4×3 in) — the "before" view |
| `best_tunnel_overlay_nolabels.png` | All profiles coloured by cluster, no legend/labels (4×3 in) — the "after" view |
| `best_tunnel_cluster_overlays.png` | Grid panel (columns ≈ √n) with the individual channel profiles per cluster only — no mean/std |
| `best_tunnel_cluster_means.png` | Grid panel (columns ≈ √n) with individual profiles + mean ± std for every cluster |
| `cluster_<N>.png` | One figure per cluster — individual profiles + mean ± std, with member list in the title |
| `cluster_<N>_overlay.png` | One figure per cluster — individual channel profiles only, no mean/std |
| `cluster_assignments.csv` | Table of `stem, tunnel, cluster, length` — the cluster label assigned to each structure, for downstream joins (see [Affinity by cluster](#affinity-by-cluster)) |

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--results-dir PATH` | `results/` | Root results folder |
| `--fpocket-dir PATH` | — | fpocket outputs folder; only PASS structures are clustered |
| `--tunnel INT` | — | Override: force this tunnel ID for all structures instead of reading `tunnel_selection.csv` |
| `--n-points INT` | `100` | Points to resample each profile to |
| `--n-clusters INT` | auto | Number of clusters to colour (default: `max(2, n // 3)`) |
| `--region LO HI` | — | Cluster using **only** the normalized-distance window `[LO, HI]` (0 = heme, 1 = exit); repeatable to combine windows; see [Region-restricted clustering](#region-restricted-clustering) |
| `--out PATH` | `results/clustering/` | Output directory |
| `--ylim FLOAT` | `5` | Y-axis half-range for overlay and mean plots |

#### Region-restricted clustering

By default the whole channel profile drives the clustering. Passing `--region LO HI` restricts the
clustering to the resampled points whose normalized position falls within `[LO, HI]` — useful when a
particular stretch of the channel (e.g. the mouth) is the biologically discriminating region. The flag
is **repeatable**: give it several times to cluster on the **union** of the windows (for example the
simple spread window *and* the smoothness window from `find_reference_region.py`). When `--region` is
set:

- the pinned entry/exit/bottleneck features are **dropped**, so the grouping depends *solely* on the
  chosen window(s);
- every profile plot still shows the **full** channel, with each clustering window **shaded grey**;
- the dendrogram title is annotated with `[regions …]`.

```bash
# single window
python3 python_scripts/cluster_tunnels.py --region 0.69 1.00 --out results/clustering_region

# two windows at once — clusters on their union
python3 python_scripts/cluster_tunnels.py \
  --region 0.08 0.24 --region 0.69 1.00 --out results/clustering_region
```

Because a region run overwrites the standard output files (including `cluster_assignments.csv`, which
the downstream scripts read), send it to a separate folder with `--out` if you want to keep it
alongside the full-profile run.

Use [`find_reference_region.py`](#finding-a-region-to-cluster-on) to choose the window(s).

---

## Finding a region to cluster on

`find_reference_region.py` uses a small set of experimentally characterized **reference structures**
(default `5OXU_heme`, `5FUJ_heme`, `9HE9_heme`) to locate a good normalized-distance region, which can
then be fed to `cluster_tunnels.py --region LO HI` to focus the clustering on a chosen stretch instead
of the whole channel.

It reuses the exact loading, normalization, and resampling of `cluster_tunnels.py`, so the profiles
are identical to those used for clustering. Because there are only a handful of references, the
analysis is deliberately **descriptive** — no ANOVA/PCA/supervised models, which are meaningless at
n ≈ 3. Two views are produced:

**Simple view (`reference_region_simple.png`)** — the region where the profiles are most **far apart**,
by per-position spread (`std` or `range`). Easy to explain, but "far apart" is not the same as
"informative".

**Smoothness view (`reference_region.png`)** — the region where the **relationships between the lines
are most coherent**. The idea: a region is low-information when the lines independently jitter up and
down (their differences are noisy), and high-information when they move together so their differences
evolve smoothly — regardless of how far apart they are (magnitude is ignored). For every pair of lines
we take the signed difference, split it into a Gaussian-smoothed trend plus a residual "wiggle", and
score `smoothness = 1 − normalized local roughness of that residual`. A flat-but-separated region
(constant difference) scores high; a region where the lines move independently scores low.

The window is taken around the peak of the chosen score (all positions ≥ `--window-frac` of the peak).
The **pairwise absolute-difference** panel is shown in both as a check that the window separates *all*
references, not just one outlier.

**Run:**

```bash
python3 python_scripts/find_reference_region.py
# then feed the reported window into region-restricted clustering:
python3 python_scripts/cluster_tunnels.py --region <LO> <HI> --out results/clustering_region
```

**Outputs** (to `results/clustering/` by default):

| File | Description |
|------|-------------|
| `reference_region_simple.png` | Simple view: reference profiles overlaid (spread window shaded), the **spread-only** curve (peak + window), and pairwise differences |
| `reference_region.png` | Smoothness view: reference profiles overlaid (smoothness window shaded); the smoothness score with spread shown for context (peak + window); and the pairwise-difference curves |
| `reference_region_scores.csv` | Per-position `spread`, `roughness`, `smoothness`, and each pairwise `|ΔR|` |

Both figures are always produced; the run prints a ready-to-use `--region` command for each window.

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--results-dir PATH` | `results/` | Root results folder |
| `--references LIST` | `5OXU_heme,5FUJ_heme,9HE9_heme` | Comma-separated reference stems |
| `--tunnel INT` | best kept | Force a tunnel ID for all references (default: best kept from `tunnel_selection.csv`) |
| `--n-points INT` | `100` | Points to resample each profile to (must match the clustering run) |
| `--spread-metric std\|range` | `std` | Magnitude-of-separation measure for the simple view |
| `--smooth-sigma FLOAT` | `3.0` | Gaussian sigma (in points) for the trend/wiggle split used by the smoothness score |
| `--window-frac FLOAT` | `0.5` | Window = contiguous region around the peak where the score ≥ this fraction of the peak |
| `--ylim FLOAT` | `5` | Y-axis half-range for the overlay panel |
| `--out PATH` | `results/clustering/` | Output directory |

> **Confirm on the pairwise panel.** Whatever window you pick, check the pairwise-difference panel that
> *all* pairwise curves are elevated in it, so it separates every reference rather than just one.

---

## Affinity by cluster

`cluster_affinity_violin.py` links docking/binding predictions to the tunnel-shape clusters and
draws violin plots of predicted affinity per cluster. It reads a [Boltz-2](https://github.com/jwohlwend/boltz)
summary table (one row per accession) and joins it to the `cluster_assignments.csv` written by
`cluster_tunnels.py`. Run `cluster_tunnels.py` **first** so that file exists.

**Pipeline:**
1. Read the summary CSV (`--summary-csv`).
2. Keep only rows whose confidence column (`--iptm-col`, default `ligand_iptm`) is **≥** the cutoff
   (`--iptm-cutoff`, default `0.8`). `ligand_iptm` is the overall complex confidence and is used
   purely as a quality gate.
3. Join on accession (first column, `--accession-col`) to the cluster assignments — only accessions
   present in **both** tables are kept.
4. Group the affinity column (`--affinity-col`, default `affinity_pred_value`) by cluster and draw
   one violin per cluster, overlaying every individual candidate as a point.

Clusters with only one member show just the point (a violin needs ≥ 2 values for its density
estimate). Each x-axis label shows the cluster's member count.

**Typical workflow:**

```bash
# Step 1 — produce cluster_assignments.csv
python3 python_scripts/cluster_tunnels.py \
  --fpocket-dir ../fpocket_pipeline/outputs

# Step 2 — violin plots of affinity across those clusters
python3 python_scripts/cluster_affinity_violin.py \
  --summary-csv boltz2_summary/boltz_nnbt_summary.csv
```

**Outputs** (to `boltz2_summary/` by default):

| File | Description |
|------|-------------|
| `affinity_violin_by_cluster.png` | Violin of `affinity_pred_value` per cluster, individual candidates overlaid as points, colours matching the tunnel clusters |
| `affinity_by_cluster.csv` | Every joined candidate (`accession, cluster, ligand_iptm, affinity_pred_value`) sorted by affinity ascending — a ranked shortlist |

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--summary-csv PATH` | *(required)* | Boltz summary CSV, one row per accession |
| `--clusters-csv PATH` | `results/clustering/cluster_assignments.csv` | Cluster labels to join against |
| `--accession-col NAME` | `accession` | Column holding the accession / structure stem (join key) |
| `--iptm-col NAME` | `ligand_iptm` | Confidence column used for filtering |
| `--iptm-cutoff FLOAT` | `0.8` | Keep rows with `iptm-col` ≥ this value |
| `--affinity-col NAME` | `affinity_pred_value` | Affinity column plotted in the violins |
| `--out PATH` | `boltz2_summary/` | Output directory |

> **Note on affinity direction:** the ranked CSV sorts ascending; confirm whether a lower
> `affinity_pred_value` corresponds to a stronger binder in your Boltz setup before reading the
> shortlist top-down.

---

## Tunnel cluster vs CUPP clade

`cluster_clade_heatmap.py` compares the **tunnel-shape clusters** (from `cluster_tunnels.py`) against
an independent **CUPP sequence clade** labelling, to ask whether structures that share a tunnel shape
also share a sequence clade. It joins the two on accession and draws a contingency heatmap. Run
`cluster_tunnels.py` **first** so `cluster_assignments.csv` exists.

**Inputs:**
- `--clusters-csv` (default `results/clustering/cluster_assignments.csv`) — `stem, cluster` from `cluster_tunnels.py`.
- `--cupp-file` (default `cupp_summary/cupp_clades.txt`) — tab-separated, one line per accession. The
  **accession** is the first field and the **clade** is the **last 3 characters** of the line (e.g. `2.2`).

Only accessions present in **both** tables are counted; any clustered structure missing from the CUPP
file is reported by name (e.g. crystal structures named by PDB ID rather than a UniProt accession).

**Run:**

```bash
python3 python_scripts/cluster_clade_heatmap.py
```

**Outputs** (to `cupp_summary/` by default):

| File | Description |
|------|-------------|
| `cluster_vs_clade_heatmap.png` | Heatmap — rows = tunnel clusters, columns = CUPP clades, each cell = number of structures with that combination (annotated with the count) |
| `cluster_vs_clade.csv` | The same contingency table with row and column totals |

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--clusters-csv PATH` | `results/clustering/cluster_assignments.csv` | Tunnel cluster labels |
| `--cupp-file PATH` | `cupp_summary/cupp_clades.txt` | CUPP clade file (accession + clade as last 3 chars) |
| `--out PATH` | `cupp_summary/` | Output directory |
| `--cmap NAME` | `Blues` | Matplotlib colormap for the heatmap |

---

## Output files (CAVER)

After a run the key files inside `out/` (or `results/<stem>/out/`) are:

| File | Description |
|------|-------------|
| `analysis/tunnel_profiles.csv` | R, X, Y, Z, distance values along each tunnel |
| `analysis/tunnel_characteristics.csv` | Bottleneck radius, throughput, length, curvature per tunnel |
| `analysis/ids_and_clusters.csv` | Mapping of tunnel IDs to clusters |
| `summary.txt` | Human-readable summary |
| `data/clusters_timeless/` | Per-cluster tunnel geometry as PDB files |
| `pymol/` | PyMOL visualisation scripts |
| `vmd/` | VMD visualisation scripts |
