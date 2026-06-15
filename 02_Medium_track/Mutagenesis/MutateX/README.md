# MutateX — In Silico Saturation Mutagenesis

# Obs! Not able to include HEME GROUP - so dont use for UPOs!

MutateX is an automated pipeline for in-silico saturation mutagenesis of protein structures using FoldX. It computes ΔΔG values for all possible single-point mutations across a protein structure, helping identify stabilizing or destabilizing mutations.

**Source:** Cancer Structural Biology Group, Danish Cancer Society Research Center / Technical University of Denmark.
**Citation:** Tiberti et al., _Brief Bioinform._ 2022 Mar 22;bbac074. doi: 10.1093/bib/bbac074

---

## Directory structure

```
MutateX/
├── mutatex/                  # MutateX source code (submodule)
├── <PROTEIN_ID>/             # One folder per protein, e.g. 6EKZ/ or A0A0C9VQZ5/
│   ├── <protein>_clean.pdb           # Cleaned input structure
│   ├── mutation_list.txt             # Amino acids to mutate TO (one-letter codes)
│   ├── position_list.txt             # (optional) Which residues to mutate
│   ├── repair_runfile_template.txt   # FoldX RepairPDB settings
│   ├── mutate_runfile_template.txt   # FoldX BuildModel settings
│   ├── rotabase.txt                  # FoldX rotamer library (copy from shared install)
│   └── run_<PROTEIN_ID>.sh           # Job script (submit to HPC)
└── README.md
```

---

## Input files

You need to prepare the following files inside your protein folder before running.

### 1. `<protein>_clean.pdb`

The input protein structure in PDB format. It must be cleaned before use:

- Remove HETATM records (ligands, waters) unless needed
  - remove solvents
  - remove (hetatm and not resn HEM)
- Keep only the chain(s) of interest
- Ensure no missing residues that FoldX cannot handle

**How to get it:** Download from RCSB PDB or use an AlphaFold model. Clean with PDB tools or manually.

### 2. `mutation_list.txt`

A plain text file listing the **target amino acids** to mutate each position _to_, one per line using single-letter codes. MutateX will try every residue in `position_list.txt` (or the whole protein if no position list is given) and mutate it to each amino acid listed here.

Example (all 20 amino acids for full saturation mutagenesis):

```
A
C
D
E
F
G
H
I
K
L
M
N
P
Q
R
S
T
V
W
Y
```

For a quick test, you can use a subset (e.g. just `G`, `A`, `V`).

### 3. `position_list.txt` _(optional but recommended)_

A plain text file specifying **which residues to mutate**. Without this file, MutateX mutates every residue in the entire protein — which is slow and often unnecessary.

Each line is one residue in the format `<AminoAcid><Chain><ResidueNumber>`, e.g.:

```
FA69
DA70
LA145
NB32
```

- The **one-letter amino acid code** of the wild-type residue comes first (e.g. `F` for Phenylalanine, `D` for Aspartate).
- Followed by the **chain letter** (e.g. `A`, `B`).
- Followed immediately by the **residue number**.
- To find the right amino acid code, check the PDB file: `grep "^ATOM" protein.pdb | awk '$6==69 && $3=="CA"' | awk '{print $4}'` then convert three-letter to one-letter code.
- For multimeric proteins where you want to mutate symmetrically linked positions, join them with `_` on a single line: `FA100_FA200`
- Residue numbers must match exactly what is in the PDB file — check with a viewer like PyMOL or by inspecting the PDB file directly.

To use the position list, add `-q position_list.txt` to your run script.

### 4. `repair_runfile_template.txt`

FoldX configuration for the `RepairPDB` step, which energy-minimizes the input structure before mutagenesis. Copy and adapt:

```
command=RepairPDB
pdb=$PDBS$
temperature=298
water=-CRYSTAL
complexWithDNA=true
```

- `$PDBS$` is a placeholder — MutateX fills it in automatically.
- Set `complexWithDNA=false` if your protein does not interact with DNA.
- `temperature=298` is 25°C (standard).

### 5. `mutate_runfile_template.txt`

FoldX configuration for the `BuildModel` step, which introduces each mutation. Copy and adapt:

```
command=BuildModel
pdb=$PDBS$
mutant-file=individual_list.txt
water=-CRYSTAL
numberOfRuns=$NRUNS$
complexWithDNA=true
```

- `$PDBS$` and `$NRUNS$` are placeholders filled by MutateX.
- `numberOfRuns` controls how many FoldX replicates are averaged (typically 5).
- Set `complexWithDNA=false` if not relevant.

### 6. `rotabase.txt`

The FoldX rotamer library. If FoldX is installed in a shared location on the HPC, you can symlink it rather than copying:

```bash
ln -s /shared/path/to/foldx/rotabase.txt rotabase.txt
```

Or copy it:

```bash
cp /shared/path/to/foldx/rotabase.txt .
```

### 7. `run_<PROTEIN_ID>.sh`

The shell script that launches MutateX — also serves as the HPC job submission script (see below).

---

## Installation on HPC - Only for Kristine right now.

If MutateX and FoldX are already installed in a **shared location on the HPC**, individual users do not need to install anything themselves — they only need to:

1. Activate the shared conda environment (or load the module), and
2. Prepare their protein folder with the input files listed above.

Ask your HPC admin for the exact paths. They will look something like:

```
Shared mutatex env:  /shared/envs/mutatex-env
FoldX binary:        /shared/software/foldx/foldx
rotabase.txt:        /shared/software/foldx/rotabase.txt
```

### Setting up the shared installation (do this once, as admin or designated person)

```bash
# 1. Create a shared conda environment
conda create -p /shared/envs/mutatex-env "python>=3.10"
conda activate /shared/envs/mutatex-env
conda install -c conda-forge setuptools biopython numpy matplotlib scipy six openpyxl pandas pyyaml
pip install logomaker "adjustText>=0.8" setuptools

# 2. Install mutatex into the shared environment
cd /path/to/MutateX/mutatex
python setup.py install

# 3. Make FoldX available (after downloading from https://foldxsuite.crg.eu)
chmod +x foldx_*
# Place foldx binary and rotabase.txt in e.g. /shared/software/foldx/
```

Optionally, set shared environment variables in the shared env's activation script so users don't have to specify paths manually:

```bash
export FOLDX_BINARY=/shared/software/foldx/foldx
export FOLDX_ROTABASE=/shared/software/foldx
```

---

## Running on HPC

MutateX is computationally intensive and **must be run on HPC** for full saturation mutagenesis.

### 1. Transfer your protein folder to HPC

```bash
scp -r <PROTEIN_ID>/ <username>@<hpc-address>:/path/to/workdir/
```

### 2. Write a job script (LSF/bsub)

Create `run_<PROTEIN_ID>.sh` inside your protein folder:

```bash
#!/bin/bash
#BSUB -J mutatex_<PROTEIN_ID>
#BSUB -n 4
#BSUB -R "rusage[mem=16GB]"
#BSUB -W 24:00
#BSUB -o mutatex_%J.log
#BSUB -e mutatex_%J.err

# Activate the shared mutatex environment
conda activate /shared/envs/mutatex-env

cd /path/to/workdir/<PROTEIN_ID>

mutatex <protein>_clean.pdb \
    -x /shared/software/foldx/foldx \
    -f suite5 \
    -p 4 \
    -m mutation_list.txt \
    -q position_list.txt \
    -R repair_runfile_template.txt \
    -M mutate_runfile_template.txt \
    -b /shared/software/foldx/rotabase.txt
```

Remove `-q position_list.txt` if you want to mutate all residues.

### 3. Submit the job

```bash
bsub < run_<PROTEIN_ID>.sh
```

### 4. Monitor your job

```bash
bjobs                 # list your running/pending jobs
bjobs -l <jobID>      # detailed info on a specific job
bpeek <jobID>         # peek at live stdout of a running job
bkill <jobID>         # cancel a job
```

**Note on resources:** Adjust `-n` (cores) to match the `-p` argument in the mutatex command. Adjust `-W` (walltime) based on protein size and how many positions/mutations you are running — a ~500 residue protein with all 20 amino acids and no position filter can take 24–48 hours.

### 5. Key MutateX arguments

| Flag         | Description                                                 |
| ------------ | ----------------------------------------------------------- |
| (positional) | Input PDB file                                              |
| `-x`         | Path to FoldX binary                                        |
| `-f`         | FoldX version: `suite4` or `suite5`                         |
| `-p`         | Number of parallel processes                                |
| `-m`         | `mutation_list.txt` — which amino acids to mutate to        |
| `-q`         | `position_list.txt` — which residues to mutate (omit = all) |
| `-R`         | Repair runfile template                                     |
| `-M`         | Mutate runfile template                                     |
| `-b`         | Path to `rotabase.txt`                                      |

---

## Post-run: cleaning up output

The `mutations/` directory produced by MutateX contains many intermediate PDB files and can be very large. After the run completes, compress it:

```bash
cd mutations
find . -name runfile.txt | xargs rm
find . -name individual_list.txt | xargs rm
find . -name rotabase.txt | xargs rm
find . -name '*.pdb' | xargs rm
cd ..
tar cvzf mutations.tar.gz mutations
rm -rf mutations
```

---

## Output

MutateX produces ΔΔG values (kcal/mol) for each mutation. Negative ΔΔG = stabilizing, positive = destabilizing. Results can be visualized with the bundled tools:

```bash
ddg2heatmap   # heatmap of ΔΔG values across the sequence
ddg2summary   # summary statistics
ddg2excel     # export to Excel
ddg2logo      # sequence logo weighted by ΔΔG
```

---

## Existing protein runs

| Protein | PDB/ID     | Folder        |
| ------- | ---------- | ------------- |
| Laccase | 6EKZ       | `6EKZ/`       |
| Laccase | A0A0C9VQZ5 | `A0A0C9VQZ5/` |

# Post-processing results

## Making a heatmap

From then same protein folder run:

'''bash
ddg2heatmap -p <input-enzyme>.pdb -d results/mutation_ddgs/final_averages -l mutation_list.txt -q position_list.txt
'''

Will output a heatmap of the change in Gibbs free energy of each mutation, as a pdf.
