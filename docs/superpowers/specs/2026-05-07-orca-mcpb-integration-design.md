# ORCA/MCPB.py Integration Fix — Design Spec
Date: 2026-05-07

## Problem

`pada1_md_pipeline.py` supports ORCA as a QM engine (stage 3), but the current
implementation is broken end-to-end:

- Stage 2 (`stage2_mcpb_step1`) always generates Gaussian `.com` input files and
  never tells MCPB.py to use ORCA, even when `qm_engine = "orca"`.
- Stage 3 (`_stage3_orca`) passes those `.com` files directly to ORCA, which
  cannot read Gaussian format.
- Stage 4 (`stage4_mcpb_steps234`) calls MCPB.py steps 2–4 to parse QM output;
  if stage 3 somehow ran, MCPB.py would fail because it expects Gaussian log files.
- Dead GAMESS-patching code in stage 2 patches Gaussian keywords into `.com` files
  without converting the format — never worked.

## Solution: MCPB.py Native ORCA Support (Approach A)

AmberTools ≥22 (available via conda) supports `quantum_soft orca` in MCPB.py
input. This causes MCPB.py to:
- Generate ORCA `.inp` files in step 1 (not Gaussian `.com`)
- Parse ORCA `.out` files in steps 2–4 (force constants, ESP charges)

This is the only approach that fixes the full pipeline without fragile format
conversion.

## Changes Required

### 1. `DEFAULT_CONFIG` — add ORCA-specific settings

```python
"orca_exe":    "/root/orca_6_1_1/orca",  # already in pada1_config.yaml
"orca_nproc":  4,     # number of parallel cores for ORCA
"orca_mem_mb": 4000,  # memory per core in MB (%maxcore in ORCA)
```

Remove `"gamess_exe"` and `"gamess_version"` from DEFAULT_CONFIG (GAMESS never
worked and is not being compiled).

### 2. `stage2_mcpb_step1` — ORCA-aware MCPB.py input

When `qm_engine == "orca"`, add to the MCPB.py `.in` file:
```
quantum_soft orca
orca_mem_mb  <cfg["orca_mem_mb"]>
```

Also when `qm_engine == "orca"`, expect `.inp` output files instead of `.com`:
- `{pdb_id}_heme_small_opt.inp`
- `{pdb_id}_heme_small_fc.inp`
- `{pdb_id}_heme_large_mk.inp`

Remove the GAMESS patching block (lines 386–393). Remove the comment claiming
MCPB.py always generates Gaussian inputs.

For Gaussian, keep patching `.com` files as before (nproc, mem, charge, mult).
For ORCA, no patching needed — MCPB.py writes the correct ORCA format directly.

### 3. `_stage3_orca` — use `.inp` paths from stage 2

Currently hardcodes path lookup from `paths["small_opt"]` etc., which point to
`.com` files. Update to use `.inp` filenames and output `.out` files:

```python
jobs = [
    (paths["small_opt"], f"{pdb_id}_heme_small_opt.out", "..."),
    ...
]
```

`small_opt`, `small_fc`, `large_mk` in the paths dict will now carry `.inp`
paths when ORCA is used (set by stage 2).

The run command is already correct: `orca {inp_file} > {outname}`.

### 4. `_recover_paths` — engine-aware path recovery

When recovering mid-pipeline, infer file extension from `cfg["qm_engine"]`:
- ORCA: look for `.inp` files
- Gaussian: look for `.com` files

### 5. Config files — remove stale GAMESS keys, confirm ORCA settings

`pada1_config.yaml` and `pada1_test_config.yaml`: ensure they have:
```yaml
qm_engine:    "orca"
orca_exe:     "/root/orca_6_1_1/orca"
orca_nproc:   4
orca_mem_mb:  4000
```

## Out of Scope

- GAMESS support (source not compiled; not needed)
- Installing AmberTools/MCPB.py (separate step: `conda install -c conda-forge ambertools`)
- OpenMPI setup for ORCA parallelism (single-core ORCA works without MPI)

## Testing

After implementation, verify with:
```bash
python pada1_md_pipeline.py --config pada1_test_config.yaml --stage 2
```
This runs only up to stage 2 and checks that `.inp` files are generated.
Full stage 3 test requires AmberTools conda env active.
