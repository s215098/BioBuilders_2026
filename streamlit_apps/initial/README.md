# 🧬 BioBuilder Suite

A multi-page Streamlit app for bioinformatics workflows including phylogenetic clustering, molecular docking, enzyme search, and more.

---

## Pages

| Page | Description | Key dependencies |
|------|-------------|-----------------|
| `phylo_clustering_app` | MSA + ML phylogenetic tree | MAFFT, FastTree, Biopython |
| `docking_app` | Molecular docking | AutoDock Vina, Boost |
| `enzyme_search_app` | Enzyme search | requests, pandas |
| `substrate_search_app` | Substrate search | requests, pandas |
| `biobuilder_info_app` | Project info | — |
| `our_project_app` | Project overview | — |

---

## Setup

### 1. Clone the repository

```bash
git clone <your-repo-url>
cd streamlit_apps/initial
```

### 2. Create the conda environment

This installs all conda-managed packages including MAFFT, PhyML, Boost, SWIG, and Python itself.

```bash
conda env create -f environment.yml
conda activate phylo_env
```

> If you already have the environment and want to update it after changes to `environment.yml`:
> ```bash
> conda env update -f environment.yml --prune
> ```

### 3. Install system binaries

**macOS (Homebrew):**
```bash
brew install fasttree
```

**Linux (apt):**
```bash
sudo apt install fasttree
```

> MAFFT and PhyML are already handled by conda in step 2.
> FastTree must be installed separately as it isn't available in the conda environment on all platforms.

### 4. Install Python packages

```bash
pip install -r requirements.txt
```

### 5. Run the app

```bash
streamlit run initial_app.py
```

---

## Project structure

```
streamlit_apps/initial/
├── initial_app.py              ← main entry point
├── requirements.txt            ← pip packages (all subapps)
├── packages.txt                ← system binaries for Streamlit Cloud
├── environment.yml             ← full conda environment snapshot
├── README.md                   ← this file
├── biobuilder_info_app/
├── docking_app/
├── enzyme_search_app/
├── our_project_app/
├── phylo_clustering_app/
└── substrate_search_app/
```

---

## Dependency overview

| Type | File | Used by |
|------|------|---------|
| Python packages | `requirements.txt` | pip, Streamlit Cloud |
| System binaries | `packages.txt` | Streamlit Cloud (apt) |
| Full conda env | `environment.yml` | local development |

---

## Updating the environment

After installing new packages, update the relevant files:

```bash
# After pip install something
pip freeze > requirements.txt

# After conda install something
conda env export > environment.yml

# Remember to remove the last 'prefix:' line from environment.yml —
# it contains your local path and will cause errors on other machines.
```

---

## Streamlit Cloud deployment

Streamlit Cloud automatically reads:
- `requirements.txt` → pip installs all Python packages
- `packages.txt` → apt installs all system binaries

Both files must be in the project root.

---

## Troubleshooting

**FastTree not found:**
```bash
brew install fasttree       # macOS
sudo apt install fasttree   # Linux
which FastTree              # verify it's on PATH
```

**MAFFT not found:**
```bash
conda install -c bioconda mafft
which mafft
```

**Vina / Boost install fails:**
```bash
conda install -c conda-forge boost
pip install vina==1.2.7
```

**Conda environment conflicts:**
```bash
conda env remove -n phylo_env
conda env create -f environment.yml
```