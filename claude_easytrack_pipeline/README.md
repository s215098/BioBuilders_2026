# EasyTrack Drylab Pipeline
### BioBuilders iGEM 2026

A command-line pipeline for discovering and evaluating enzyme candidates for a given substrate. Starting from a keyword search, it builds a curated sequence library, constructs a phylogenetic tree to select diverse representatives, and runs AutoDock Vina molecular docking to predict how each candidate binds the substrate of interest.

---

## Pipeline Overview

```
  ┌─────────────────────────────────────────────────────────────┐
  │  Config file  (config/your_enzyme_ligand.yml)               │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 1 │ Sequence Retrieval                                │
  │         │ NCBI protein search  OR  UniProt API              │
  │         │ → sequences/raw_sequences.fasta                   │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 2 │ MSA + Phylogenetic Tree                           │
  │         │ MAFFT alignment  →  FastTree (ML tree)            │
  │         │ → phylogeny/aligned.fasta                         │
  │         │ → phylogeny/tree.nwk  +  tree.png                 │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 3 │ Select Representatives                            │
  │         │ Cluster tree by branch-length distance            │
  │         │ Pick one diverse representative per cluster       │
  │         │ → phylogeny/clusters.tsv  +  selections.txt       │
  │         │ → sequences/representatives.fasta                 │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 4 │ Receptor Preparation                              │
  │         │ Download PDB from RCSB  (or use local file)       │
  │         │ pdbfixer: clean + add missing atoms + protonate   │
  │         │ OpenBabel: convert to PDBQT                       │
  │         │ → receptor/<ID>_clean_h.pdbqt                     │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 5 │ Ligand                                            │
  │         │ Download 3D SDF from PubChem by name or CID      │
  │         │ → ligand/<name>.sdf                               │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 6 │ Molecular Docking (AutoDock Vina)                 │
  │         │ Blind docking (no prior knowledge needed)         │
  │         │ OR targeted (known active site coordinates)       │
  │         │ → docking/<run>/poses.pdbqt  +  poses.pdb         │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 6b│ Boltz-2 + Targeted Vina  [OPTIONAL]               │
  │         │ Boltz-2 predicts complex from sequence + SMILES   │
  │         │ → extracts AI-predicted pocket center             │
  │         │ → runs targeted Vina using that pocket            │
  │         │ → boltz/<run>/boltz_summary.json                  │
  │         │ → boltz/<run>/final_poses.pdbqt                   │
  │         │ Enable with:  boltz.enabled: true  in config      │
  └──────────────────────┬──────────────────────────────────────┘
                         │
  ┌──────────────────────▼──────────────────────────────────────┐
  │  Step 7 │ Parse Results                                     │
  │         │ Extract binding affinities from all poses         │
  │         │ → docking/summary.csv  +  best_poses.txt          │
  └─────────────────────────────────────────────────────────────┘
```

---

## Installation

> **Python 3.12 is required.** The `vina` Python package has pre-built wheels for 3.12 but not for 3.10/3.11. Building from source requires Boost headers to be explicitly set.

### Step-by-step

```bash
# 1. Create and activate the conda environment
#    (installs MAFFT, FastTree, pdbfixer, OpenBabel, and build deps)
conda env create -f environment.yml
conda activate easytrack

# 2. Export Boost paths so pip can build/find vina
#    (these must be set BEFORE pip install vina)
export BOOST_INCLUDE=$CONDA_PREFIX/include
export BOOST_LIB=$CONDA_PREFIX/lib

# 3. Install the Python packages that need pip
pip install vina==1.2.7 meeko==0.7.1 rdkit==2025.9.3
```

> **If you already created the env and vina fails to import:**
> ```bash
> conda activate easytrack
> conda install -c conda-forge swig boost-cpp -y
> export BOOST_INCLUDE=$CONDA_PREFIX/include
> export BOOST_LIB=$CONDA_PREFIX/lib
> pip install vina==1.2.7
> ```

### Verify everything works

```bash
mafft --version
FastTree -help 2>&1 | head -1
obabel --version

# On macOS + conda, use 'python' not 'python3' (python3 may resolve to system Python)
python -c "from vina import Vina; from rdkit import Chem; from meeko import MoleculePreparation; print('Vina + RDKit + Meeko OK')"
python -c "from pdbfixer import PDBFixer; print('pdbfixer OK')"
python -c "from Bio import Entrez, Phylo; print('Biopython OK')"
```

> **macOS tip:** Even with a conda env active, `python3` can still resolve to `/usr/bin/python3` (the system Python). Always use `python` (no 3) inside a conda env on macOS, or verify with `which python`.

---

## Quick Start — Example (8RNJ + BADGE)

```bash
conda activate easytrack
cd claude_easytrack_pipeline/

# Run the full pipeline with the provided example config
bash run_pipeline.sh config/example_8rnj_badge.yml
```

Results are written to `results/8rnj_badge/`.

---

## Configuration

All settings for a pipeline run are stored in a single YAML file. Copy the example and edit it for your enzyme and ligand:

```bash
cp config/example_8rnj_badge.yml config/my_enzyme_myligand.yml
# Edit with any text editor
nano config/my_enzyme_myligand.yml
```

Key settings to change:

| Section | Key | Description |
|---|---|---|
| `sequences` | `mode` | `ncbi`, `uniprot`, or `manual` |
| `sequences` | `query` | Search string, e.g. `"laccase AND bacteria"` |
| `sequences` | `ncbi_email` | **Required** by NCBI — use your real email |
| `selection` | `mode` | `auto` or `manual` (pause for user review) |
| `receptor` | `pdb_id` | RCSB PDB ID, e.g. `"6T0Y"` |
| `ligand` | `name` | Compound name for PubChem search |
| `docking` | `mode` | `blind` or `targeted` |
| `docking` | `center` | `[x, y, z]` for targeted docking |

---

## Running Individual Steps

Each script in `pipeline/` can be run standalone:

```bash
# On macOS + conda, use 'python' not 'python3'
python pipeline/01_fetch_sequences.py --config config/example_8rnj_badge.yml
python pipeline/02_build_phylogeny.py --config config/example_8rnj_badge.yml
python pipeline/03_select_representatives.py --config config/example_8rnj_badge.yml
python pipeline/04_prepare_receptor.py --config config/example_8rnj_badge.yml
python pipeline/05_fetch_ligand.py --config config/example_8rnj_badge.yml
python pipeline/06_docking.py --config config/example_8rnj_badge.yml
python pipeline/07_parse_results.py --config config/example_8rnj_badge.yml
```

### Resume from a specific step

```bash
# Restart from step 4 onwards (skips sequence retrieval + phylogeny)
bash run_pipeline.sh config/example_8rnj_badge.yml --from-step 4

# Run only docking
bash run_pipeline.sh config/example_8rnj_badge.yml --only-step 6
```

---

## Manual Representative Selection (Step 3)

Set `selection.mode: "manual"` in your config. The pipeline will:

1. Cluster the tree and generate a default `selections.txt` with one sequence per cluster
2. **Pause** and print instructions for you to review
3. Wait for you to edit `selections.txt` and re-run

```bash
# After editing phylogeny/selections.txt:
bash run_pipeline.sh config/my_config.yml --from-step 3 --apply-selections
```

The file `phylogeny/clusters.tsv` lists all sequences and their cluster IDs, so you can pick any member from a cluster as the representative.

---

## Manual Receptor Preparation (Step 4 alternative)

If you prefer to use PyMOL for cleaning instead of the automated pdbfixer approach:

1. Open the raw PDB in PyMOL
2. Run the following commands:
```
remove chain B        # keep only chain A (adjust as needed)
remove solvent        # remove water
remove hetatm         # remove ligands and ions
h_add                 # add hydrogens
save receptor_clean_h.pdb
```
3. Convert to PDBQT with OpenBabel:
```bash
obabel -ipdb receptor_clean_h.pdb -opdbqt -O receptor_clean_h.pdbqt -xr
```
4. Place the PDBQT in `results/<your_run>/receptor/` and continue from step 5:
```bash
bash run_pipeline.sh config/my_config.yml --from-step 5
```

---

## Boltz-2 + Targeted Vina (Step 6b)

This optional step combines two complementary approaches for a richer result:

**Stage 1 — Boltz-2 AI prediction:** Given only the protein sequence and ligand SMILES, Boltz-2 predicts the full protein-ligand complex structure and estimates binding affinity as log10(IC50). No crystal structure needed.

**Stage 2 — Targeted Vina refinement:** The predicted ligand position is extracted as a pocket center, which is then used to run a focused AutoDock Vina search — more accurate than blind docking because the search space is informed by Boltz.

### Setup

```bash
pip install boltz    # ~2 GB model download on first run
```

### Config

In your config file, add:
```yaml
boltz:
  enabled: true
  sequence: "MSTAGKVIKCKAAVLWEE..."   # full amino acid sequence
  smiles: "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1"  # BADGE
  confidence_threshold: 0.60
```

Then run step 6b directly:
```bash
python pipeline/06b_boltz_vina.py --config config/my_config.yml
```

Or it runs automatically as part of `run_pipeline.sh` when `boltz.enabled: true`.

### Outputs

```
results/<run>/boltz/<name>/
  boltz_input.yaml        — Boltz input
  boltz_out/              — raw Boltz prediction files
  receptor_clean.pdb      — receptor with AI ligand removed
  receptor.pdbqt          — prepared for Vina
  final_poses.pdbqt       — targeted Vina poses
  boltz_summary.json      — confidence, log10(IC50), Vina affinity
```

---

## Docking Modes

Three modes are available — pick based on what you know about the binding site.

**`blind`** — searches the entire protein surface. Use when you have no prior knowledge of the binding site, or want to verify that docking finds the expected site on its own. Slowest.

**`targeted`** — focuses the search on exact Å coordinates. Use when you have a crystal structure with a co-crystallised ligand and can read the coordinates directly.

```yaml
docking:
  mode: "targeted"
  center: [18.8, 32.2, 14.9]   # read from PyMOL bottom toolbar
  box_size: [20.0, 20.0, 20.0]
```

To get coordinates in PyMOL: load the structure, click an active site atom (e.g. the heme iron), and read `[x, y, z]` from the bottom toolbar.

**`residues`** — calculates the box centre from named residues you specify. Best when the active site is described in the literature by residue name/number but you don't have exact coordinates. More portable than hardcoded coordinates because residue numbers stay consistent across homologous structures.

```yaml
docking:
  mode: "residues"
  active_site_residues: ["HEM", "CYS36", "HIS86", "ARG189"]
  box_size: [20.0, 20.0, 20.0]
```

You can mix residue names, name+number, or plain numbers:
```yaml
active_site_residues: ["HEM", "HIS86", "36"]   # all three formats work
```

Typical active site residues for your enzyme families:

| Enzyme | Key residues |
|---|---|
| UPO (heme-thiolate) | `HEM`, `CYS36`, `HIS86`, `ARG189` (adjust numbering per structure) |
| Laccase (T1 copper) | `HIS395`, `HIS447`, `CYS450`, `HIS452` (adjust per structure) |

You can also override the mode from the command line without editing the config:
```bash
python pipeline/06_docking.py --config config/my.yml --residues HEM HIS86 ARG189
python pipeline/06_docking.py --config config/my.yml --center 18.8 32.2 14.9
```

## Ligand Format — What Happens Automatically

The pipeline only needs an SDF file as input (from PubChem via step 5, or your own). The docking script handles all conversion automatically:

1. **Reads** the SDF with RDKit
2. **Adds hydrogens** (essential for accurate docking)
3. **Generates a 3D conformation** with ETKDGv3 if the SDF has no 3D coordinates
4. **Energy-minimises** with MMFF94
5. **Converts to PDBQT** with Meeko (the format Vina requires)

You will see each step confirmed in the terminal output:
```
[Ligand] Molecule: formula=C21H24O4, MW=340.4, heavy atoms=25
[Ligand] Hydrogens added (49 total atoms)
[Ligand] 3D coordinates found in SDF — using existing conformation.
[Ligand] MMFF94 energy minimisation: converged.
[Ligand] Converting to PDBQT with Meeko …
[Ligand] ✓ Ligand ready for Vina.
```

No manual conversion with PyMOL or OpenBabel is needed for the ligand.

---

## Output Structure

After a full run with `config/example_8rnj_badge.yml`:

```
results/8rnj_badge/
├── sequences/
│   ├── raw_sequences.fasta          # All downloaded sequences
│   └── representatives.fasta        # Curated diverse subset
├── phylogeny/
│   ├── aligned.fasta                # MAFFT alignment
│   ├── tree.nwk                     # Newick tree file
│   ├── tree.png                     # Tree image (open to review)
│   ├── clusters.tsv                 # All sequences with cluster IDs
│   └── selections.txt               # Editable representative list
├── receptor/
│   ├── 8RNJ_raw.pdb                 # Downloaded from RCSB
│   ├── 8RNJ_clean_h.pdb             # Cleaned + protonated
│   └── 8RNJ_clean_h.pdbqt           # Ready for Vina
├── ligand/
│   └── BADGE.sdf                    # 3D conformer from PubChem
└── docking/
    ├── 8RNJ_BADGE_blind/
    │   ├── poses.pdbqt              # All docked poses
    │   ├── poses.pdb                # Same, for PyMOL
    │   └── docking_log.txt          # Binding affinities per pose
    ├── summary.csv                  # All runs × all poses
    └── best_poses.txt               # Ranked best binding affinities
```

---

## Visualising Results in PyMOL

```bash
# Load receptor + docked poses together
pymol results/8rnj_badge/receptor/8RNJ_clean_h.pdb \
      results/8rnj_badge/docking/8RNJ_BADGE_blind/poses.pdb
```

In PyMOL:
- The receptor appears as the protein structure
- Each `poses_*.pdb` model is one docking pose
- Use `Next/Previous` in the sequence viewer to cycle through poses
- Check binding affinities in `docking/best_poses.txt`

---

## Adapting for a New Enzyme + Ligand

1. Copy the example config:
```bash
cp config/example_8rnj_badge.yml config/my_laccase_bpa.yml
```

2. Edit the key fields:
```yaml
output_dir: "results/my_laccase_bpa"
sequences:
  query: "laccase geobacillus"
receptor:
  pdb_id: "6T0Y"
ligand:
  name: "BPA"
  pubchem_cid: 6623
```

3. Run:
```bash
bash run_pipeline.sh config/my_laccase_bpa.yml
```

---

## Troubleshooting

**`mafft: command not found`**
```bash
conda install -c bioconda mafft
```

**`FastTree: command not found`**
```bash
conda install -c bioconda fasttree
# The pipeline will fall back to neighbor-joining automatically if FastTree is missing.
```

**`obabel: command not found`**
```bash
conda install -c conda-forge openbabel
```

**`pdbfixer` / `openmm` errors**
```bash
conda install -c conda-forge pdbfixer openmm
```

**`ModuleNotFoundError: No module named 'vina'`**

vina needs Boost headers visible at install time. Run:
```bash
conda install -c conda-forge swig boost-cpp -y
export BOOST_INCLUDE=$CONDA_PREFIX/include
export BOOST_LIB=$CONDA_PREFIX/lib
pip install vina==1.2.7
```
Also make sure you are using **Python 3.12** — check with `python3 --version`. If not, recreate the env: `conda env remove -n easytrack && conda env create -f environment.yml`.

**NCBI returns no results**
- Check your query string at https://www.ncbi.nlm.nih.gov/protein/
- Try fewer keywords (e.g. `"peroxygenase"` instead of `"unspecific peroxygenase AND fungi"`)
- Set `sequences.mode: "manual"` and provide your own FASTA

**No 3D conformer on PubChem**
- The script downloads a 2D SDF and RDKit generates the 3D conformer during docking
- Alternatively, search at https://pubchem.ncbi.nlm.nih.gov/ and download manually, then set `ligand.local_sdf`
