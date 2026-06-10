# Boltz-2 structure prediction — generic template-driven batch

Self-contained folder for predicting holo structures of any protein + arbitrary cofactors / ligands using Boltz-2. **Structures only** — no docking, no ranking, no classification. Just hand it a FASTA and a YAML template, get one PDB per sequence + a metrics CSV.

```
boltz_structures/
├── README.md
├── boltz_template_batch.py     # the batch driver
└── templates/
    ├── heme.yaml               # protein + heme b (affinity on heme)
    ├── heme_nnbt.yaml          # protein + heme + NNBT (affinity on NNBT)
    └── 4cu.yaml                # protein + 4 Cu(II) ions (no substrate)
```

---

## 1. Environment

```bash
conda create -n docking python=3.12 -y
conda activate docking
pip install boltz
```

For CUDA-13 / torch ≥ 2.12 envs the script passes `--no_kernels` to Boltz so it skips the `cuequivariance-torch` kernel path (which only ships CUDA-12 wheels). Functionally equivalent, just slightly slower.

Confirm GPU:
```bash
python -c "import torch; print('cuda ok:', torch.cuda.is_available())"
```

---

## 2. Run
Rememeber to run SignalP on your fasta file, https://services.healthtech.dtu.dk/services/SignalP-6.0/ to remove the signal peptides.
```bash
python boltz_template_batch.py \
    -f your_sequences.fasta \
    -t templates/4cu.yaml \
    -o outputs_4cu \
    --resume
```

Background-friendly:
```bash
nohup python boltz_template_batch.py -f your_sequences.fasta -t templates/4cu.yaml -o outputs_4cu --resume > boltz.log 2>&1 &
```

Required flags:
- `-f / --fasta` — input FASTA. Accession = first whitespace token after `>`.
- `-t / --template` — YAML template (see templates/ for examples, or write your own).
- `-o / --out_dir` — where to drop outputs.

Useful optional flags:
- `--resume` — skip accessions whose structure already exists; re-collect metrics for them.
- `--pdb-dir` — override the collected-PDB folder (default `<out_dir>/structures_pdb`).
- `--csv` — override the summary CSV path.
- `--tag` — override the log prefix (default derived from template name, e.g. `BOLTZ-4CU`).

---

## 3. Writing your own template

A template is plain Boltz YAML with one mandatory placeholder, `{seq}`, where the FASTA sequence gets substituted. Optionally `{accession}` is also substituted with the FASTA accession.

Minimal protein-only template:
```yaml
sequences:
  - protein:
      id: A
      sequence: "{seq}"
```

Add cofactors as ligand blocks with CCD codes (full PDB chemical-component dictionary, e.g. `HEM` heme b, `CU` Cu²⁺, `CU1` Cu⁺, `ZN` Zn²⁺, `MG` Mg²⁺, `NAD` NAD⁺, `FAD`, `FMN`, etc.):
```yaml
  - ligand:
      id: B
      ccd: HEM
```

Add small-molecule substrates / binders with SMILES:
```yaml
  - ligand:
      id: C
      smiles: 'CC1=CC=C(C=C1)N(CC(C)O)CC(C)O'
```

Tell Boltz to predict affinity for one ligand (only one binder allowed):
```yaml
properties:
  - affinity:
      binder: C
```

The substitution uses plain `str.replace`, so YAML braces elsewhere in the file are safe. The script sanity-checks that the template has a `sequences:` block and the `{seq}` placeholder before launching.

---

## 4. Output

```
<out_dir>/
├── boltz_<template-stem>_summary.csv          one row per sequence
├── structures_pdb/
│   ├── <accession>.pdb                         predicted complex, accession-named
│   └── ...
└── <accession>/
    ├── boltz_input.yaml                        what was actually folded
    └── boltz_out/                              Boltz native outputs (confidence JSON, raw PDB)
```

CSV columns (missing values just stay blank):
```
accession, status, seq_len,
confidence_score, ptm, iptm, ligand_iptm,
protein_iptm, complex_plddt, complex_iplddt,
affinity_pred_value, affinity_probability_binary,
pdb_path, structure_path
```

- `status` is `ok` or `error: ...` (per-row failures don't stop the run).
- `confidence_score` is Boltz's umbrella confidence (0–1, higher = better).
- `affinity_*` columns are populated only if your template has a `properties: affinity:` block; otherwise blank.

---
