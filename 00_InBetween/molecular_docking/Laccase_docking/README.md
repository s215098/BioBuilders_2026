# Laccase Virtual Screening Pipeline

Virtual screening of **NNBT** and **Guaiacol** against ~76 laccase-like enzyme sequences using **AutoDock Vina** at pH 7. The T1 copper active site is located per-structure by mapping three conserved coordinating residues (His443, Cys500, His505) via BLOSUM62 sequence alignment, then centering a **22 Å search box** on the sidechain-atom centroid. Each protein×ligand pair is docked in **3 independent replicates** (seeds 1–3, exhaustiveness 24) and results are filtered by a T1 geometry check to ensure only catalytically productive poses contribute to the final ranking.

---

## Directory layout

```
project/
├── laccase_docking_pipeline.py   ← main script
├── check_env.py                  ← dependency checker
├── visible_sequences.csv         ← input: Header + Sequence columns (82 entries)
├── reference.pdb                 ← reference laccase PDB with CU HETATM (T1 Cu first)
├── itl_structures/               ← manually provided PDB files
│   ├── ItL01.pdb
│   ├── ItL02.pdb
│   └── ItL03.pdb
├── structures/                   ← auto-created: AlphaFold downloads
├── trimmed/                      ← auto-created: signal-peptide-trimmed PDBs
├── prepared/                     ← auto-created: pH-7 protonated PDBQT receptors
├── ligands/                      ← auto-created: PDBQT + SDF ligand files
└── results/                      ← auto-created: docking poses + summary CSV
```

---

## Installation

### Conda (recommended)

```bash
conda create -n docking python=3.10
conda activate docking

conda install -c conda-forge biopython rdkit numpy
conda install -c conda-forge vina          # Vina Python bindings + CLI
pip install meeko pdb2pqr                  # PDBQT prep + pH-7 protonation
conda install -c conda-forge openbabel     # last-resort fallback only
```

### pip-only

```bash
pip install biopython rdkit-pypi numpy requests meeko pdb2pqr
# Vina binary: https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

### Check everything is working

```bash
python check_env.py
```

---

## Reference PDB requirements

The reference PDB must have the **T1 Cu listed first** among HETATM CU records:

```
HETATM  NNN  CU  CU  A  701      23.347   0.050  39.904  1.00 10.66  ...
```

The script identifies the T1 site via three conserved coordinating residues hardcoded as sequence indices in `_REF_T1_RESIDUES` (His443 / Cys500 / His505 of the reference). If you swap the reference PDB, update those indices to match.

---

## Usage

### Dry run — validates structure downloads, trimming, and T1 localisation without calling Vina

```bash
python laccase_docking_pipeline.py --dry-run
```

### Full run with positive control

```bash
python laccase_docking_pipeline.py --positive-control
```

### All options

```
--csv               Input CSV (default: visible_sequences.csv)
--reference         Reference PDB with T1 Cu (default: reference.pdb)
--itl-dir           Dir with ItL*.pdb files (default: itl_structures)
--box-size          Search box full-edge in Å (default: 22)
--exhaustiveness    Vina exhaustiveness per replicate seed (default: 24)
--n-poses           Poses to save per run (default: 9)
--cpus              CPU threads for Vina (default: all available)
--dry-run           Plan only, no Vina calls
--positive-control  Dock Guaiacol into reference.pdb first as a sanity check
--struct-dir        AlphaFold cache dir (default: structures)
--trimmed-dir       Trimmed PDB dir (default: trimmed)
--prepared-dir      PDBQT receptor dir (default: prepared)
--ligand-dir        PDBQT ligand dir (default: ligands)
--results-dir       Output dir (default: results)
```

---

## Pipeline overview

1. **Structure acquisition** — AlphaFold API download for UniProt IDs; manual PDBs for ItL01/02/03. 6 IDs have no AlphaFold entry and are skipped (A0A6A6ML18, A0A4D6KPY7, A0A7J6HXP0, A0A498JA70, A0A5N6R8J3, J3M898).
2. **Signal-peptide trim** — N-terminal residues outside the reference-aligned core with pLDDT < 70 are stripped from AlphaFold models before prep.
3. **Receptor prep** — `pdb2pqr --ph 7.0` (PROPKA) assigns protonation states and His tautomers; meeko writes the PDBQT. Falls back to meeko-only → MGLTools → obabel (quality warning issued).
4. **Ligand prep** — RDKit ETKDGv3 + MMFF conformer; meeko PDBQT. Reactive atoms identified by SMARTS: guaiacol → phenolic O, NNBT → aniline N.
5. **T1 localisation** — BLOSUM62 sequence alignment maps His443/Cys500/His505 of the reference to each target; box centre = centroid of Nδ1/Nε2 (His) and Sγ (Cys) atoms.
6. **Docking** — 3 replicates per protein×ligand (seeds 1–3, exhaustiveness 24, 22 Å box). Poses from seed 1 saved to PDBQT.
7. **Geometry filter** — reactive atom distance to T1 centroid measured for each pose. `productive_pose = True` if any pose ≤ 8 Å. `ranking_score` = best affinity among productive poses (falls back to overall best if none).

---

## Outputs

| File | Description |
|------|-------------|
| `results/<ID>_<ligand>_poses.pdbqt` | Docked poses from seed 1 |
| `results/<ID>_<ligand>.seed2.pdbqt` | Poses from seed 2 |
| `results/<ID>_<ligand>.seed3.pdbqt` | Poses from seed 3 |
| `results/<ID>_<ligand>.log` | Vina affinity log (seed 1) |
| `results/docking_summary.csv` | One row per protein×ligand pair |
| `results/_positive_control_Guaiacol_poses.pdbqt` | Positive control poses (if run) |

### Summary CSV columns

| Column | Description |
|--------|-------------|
| `protein_id` | UniProt ID or ItL label |
| `ligand` | NNBT or Guaiacol |
| `best_affinity_mean` | Mean best-pose affinity across 3 seeds (kcal/mol) |
| `best_affinity_sd` | SD across seeds — high values flag noisy runs |
| `best_affinity_kcal_mol` | Best pose from seed 1 |
| `all_affinities_kcal_mol` | Semicolon-separated per-pose affinities (seed 1) |
| `t1_cu_x/y/z` | Box centre used (Å, in the receptor's coordinate frame) |
| `reactive_atom_dist_to_T1` | Distance of top pose's reactive atom to T1 centroid (Å) |
| `productive_pose` | True if any pose has reactive atom ≤ 8 Å from T1 |
| `ranking_score` | Sort key: best affinity among productive poses (most negative = best) |
| `poses_pdbqt` | Path to seed-1 poses file |

### Ranking results

```python
import pandas as pd
df = pd.read_csv("results/docking_summary.csv")
df["ranking_score"] = pd.to_numeric(df["ranking_score"], errors="coerce")
# Productive poses only, best binders first
print(df[df["productive_pose"] == True].sort_values("ranking_score").head(20))
```

---

## Visualising poses

```bash
pymol prepared/Q9SU40.pdbqt results/Q9SU40_NNBT_poses.pdbqt
```

---

## Performance

| Stage | Estimate |
|-------|----------|
| AlphaFold downloads (76 structures) | ~5 min |
| Trim + prep (76 receptors, pdb2pqr+meeko) | ~10–20 min |
| Docking (152 pairs × 3 seeds, 8 CPUs) | ~6–12 hours |

To reduce runtime: lower `--exhaustiveness` to 8 for a first-pass screen, then re-dock top hits at 24. For GPU acceleration, Gnina accepts the same PDBQT format.
