# Automated Rational Mutation Scan Pipeline

Iterative pipeline for evaluating point mutations on an enzyme–substrate pair using **Boltz-2** co-folding and affinity prediction. Claude (optional) adds structural reasoning to each step report when available.

Built by iGEM 2026 DTU BioBuilders for epoxy resin degradation, but designed to work with any enzyme and substrate.

---

## How it works

You provide a ranked list of mutations. The pipeline tests them one by one:

```
For each mutation in your list:
  1. Apply mutation to the current best sequence (string operation)
  2. Run Boltz-2  →  co-folded structure + predicted affinity + contact residues
  3. Compare Δaffinity to the previous best
     - Δ ≤ threshold  →  accept, carry forward
     - Δ > threshold  →  reject, revert to previous sequence
  4. If Claude CLI is available: add structural commentary to the step report
  5. Save step report (JSON) + running summary (CSV)
```


---

## Requirements

### Python packages

```bash
pip install boltz biopython requests pyyaml pandas numpy
```

### Boltz-2 weights

Downloaded automatically on first run (~10 GB). Or trigger manually:

```bash
boltz predict --help
```

### GPU (strongly recommended)

| Hardware | Time per prediction |
|----------|-------------------|
| CPU only | ~2–4 hours |
| RTX 3080 / 4080 | ~20–60 seconds |
| A100 / H100 | ~5–15 seconds |

No local GPU? Options:
- [Google Colab](https://colab.research.google.com) (free T4)
- [Tamarind Bio](https://app.tamarind.bio/boltz) (web interface for Boltz)
- [DTU HPC](https://docs.google.com/document/d/1JnY5QevqU_iZuEjYxGVdGpGGQCkvsjftTunuR64RGZs/edit?usp=drive_link)

### Claude CLI (optional)

If the `claude` CLI is installed and logged in, it adds structural reasoning to each step report. The pipeline runs identically without it — you just get Boltz-2 numbers only.

Install: https://claude.ai/code

---

## Quick start

### 1. Set up your config

Copy the example and edit it for your system:

```bash
cp configs/pada1_nnbt.yaml configs/my_enzyme.yaml
```

Open `configs/my_enzyme.yaml` and fill in:
- `enzyme.name` and `enzyme.sequence` (mature protein sequence, no signal peptide)
- `enzyme.protected_residues` (positions that must never be mutated)
- `substrate.name` and `substrate.smiles`
- `pocket_residues` (active site residues for Boltz-2 pocket constraint)
- `cofactor` block — remove entirely for cofactor-free enzymes
- `context` — free-text description of the reaction mechanism and key residues (only needed if using Claude)

### 2. Run a mutation scan

```bash
python3 pipeline.py \
    --config configs/my_enzyme.yaml \
    --scan_mutations A316P,F191L,A73T \
    --workdir results/my_run/
```

Mutations are tested in the order given. Each one builds on the previously accepted sequence.

### 3. Adjust the acceptance threshold (optional)

```bash
python3 pipeline.py \
    --config configs/my_enzyme.yaml \
    --scan_mutations A316P,F191L,A73T \
    --rollback_threshold 0.15 \
    --workdir results/my_run/
```

`--rollback_threshold` sets the maximum allowed Δaffinity for a mutation to be accepted. Default is `0.3` log₁₀(IC₅₀), which allows mutations that are up to ~2× worse in affinity (they may still improve geometry or selectivity). Set it lower for stricter acceptance.

---

## Output

```
results/my_run/
├── scan_summary.csv              ← one row per mutation, full trajectory
├── scan_final_sequence.txt       ← final accepted sequence in FASTA format
├── round0_WT.yaml                ← Boltz-2 input for baseline
├── round1_A316P.yaml
├── scan_step1_A316P/
│   └── step_report.json          ← Boltz-2 numbers (+ Claude commentary if available)
├── scan_step2_F191L/
│   └── step_report.json
└── boltz_out/
    └── predictions/
        ├── *.cif                 ← predicted structures (open in PyMOL or ChimeraX)
        └── affinity_*.json       ← raw Boltz-2 affinity output
```

### Reading the numbers

| Field | Meaning | What to look for |
|-------|---------|-----------------|
| `affinity_pred_value` | log₁₀(IC₅₀ in µM) | Lower = tighter binder. −1.0 ≈ 0.1 µM, +2.0 ≈ 100 µM |
| `binder_probability` | P(this is a binder) | Use for relative ranking; ~0.4 = uncertain |
| `confidence_score` | Boltz-2 structure confidence | >0.8 = reliable prediction |
| `delta_affinity` | Change vs. previous best | Negative = improvement |
| `contacts` | Residues within 4.5 Å of substrate | Compare across steps to track pose changes |

**Important caveats:**
- Differences of <0.1 in `affinity_pred_value` are likely within prediction noise.
- `binder_probability` ~0.4 means the model is uncertain — use as relative ranking, not absolute.
- Boltz-2 affinity does not capture catalytic geometry. A mutation that improves predicted affinity may still position the substrate non-productively. If your reaction has a specific geometry requirement (e.g. α-C to Fe distance for UPOs), verify it separately.

---

## Adapting to a new enzyme or substrate

Everything specific to PaDa-I and NNBT lives in `configs/pada1_nnbt.yaml`. To run on a different system:

1. **Copy the config:**
   ```bash
   cp configs/pada1_nnbt.yaml configs/my_enzyme_mysubstrate.yaml
   ```

2. **Fill in your sequence.** Use the mature protein sequence (post signal peptide). Boltz-2 numbers from position 1 — verify that residue numbers in your literature match your sequence string.

3. **Set `protected_residues`.** These positions will never be mutated. Include catalytic residues, metal-binding cysteines/histidines, disulfide pairs, etc.

4. **Set `pocket_residues`.** These bias Boltz-2 to fold the substrate into your active site. Use residues you know contact the substrate from crystal structures or literature.

5. **Remove the `cofactor` block** if your enzyme has no cofactor. Leave it in (updating `ccd` or `smiles`) if it does — common CCD codes: `HEM` (heme b), `FAD`, `NAD`, `FMN`, `ZN`.

6. **Write a `context`** (only needed for Claude). Describe: what reaction you want, what the correct binding geometry looks like, and which residues are important and why. Claude will use this to comment on each Boltz-2 result. The more specific, the better.

---

## Config reference

```yaml
enzyme:
  name: MyEnzyme              # used in log output and FASTA headers
  sequence: "MAST..."         # mature sequence, single-letter, no spaces
  protected_residues: [1, 2]  # 1-indexed positions, never mutated

substrate:
  name: MySubstrate
  smiles: "C1CCCCC1"          # SMILES string

pocket_residues: [10, 50, 100]  # optional, remove to let Boltz-2 decide freely

cofactor:                     # optional — remove for cofactor-free enzymes
  chain_id: C
  ccd: HEM                    # CCD code (preferred)
  # smiles: "..."             # or SMILES if non-standard

context: |                    # optional, required only when using Claude
  Target reaction: ...
  Key residues: ...
```

---

## Example results — PaDa-I × NNBT

Wild-type baseline (AutoDock Vina): **−6.64 kcal/mol**, α-C to Fe = 3.45 Å (borderline non-productive).

Boltz-2 scan (A316P → F191L → A73T):

| Step | Variant | Δ Affinity | Pred. affinity | P(binder) | Accepted |
|------|---------|-----------|----------------|-----------|---------|
| 0 | WT | — | 1.42 | 0.40 | baseline |
| 1 | A316P | −0.044 | 1.38 | 0.43 | yes |
| 2 | A316P+F191L | +0.043 | 1.42 | 0.38 | yes |
| 3 | A316P+F191L+A73T | +0.017 | 1.44 | 0.40 | yes |

All three accepted. The affinity differences are small (~0.06 log units total). Contact residue changes suggest real ligand repositioning around the entrance gate after F191L, partially restored after A73T.

---

## References

- Linde et al. (2022) *Front. Catal.* [doi:10.3389/fctls.2022.883263](https://doi.org/10.3389/fctls.2022.883263) — NNBT N-dealkylation, 12% WT baseline
- Ramirez-Escudero et al. (2018) — PaDa-I active site (PDB: [6EKZ](https://www.rcsb.org/structure/6EKZ))
- Molina-Espeja et al. (2014) — PaDa-I directed evolution
- [Boltz-2](https://github.com/jwohlwend/boltz) — co-folding and affinity prediction
