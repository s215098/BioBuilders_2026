# ORCA/MCPB.py Integration Fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `pada1_md_pipeline.py` so that ORCA 6.1.1 works end-to-end as the QM engine via MCPB.py's native `quantum_soft orca` support.

**Architecture:** Tell MCPB.py to generate ORCA `.inp` files directly (stage 2), run them with ORCA (stage 3), and let MCPB.py parse ORCA `.out` files natively (stage 4). No format conversion needed. Remove dead GAMESS patching code.

**Tech Stack:** Python 3, ORCA 6.1.1 at `/root/orca_6_1_1/orca`, AmberTools ≥22 (conda), pytest

---

## File Map

| File | Change |
|------|--------|
| `MD/pada1_md_pipeline.py` | Modify DEFAULT_CONFIG, stage2, _recover_paths |
| `MD/pada1_config.yaml` | Add orca_nproc / orca_mem_mb, fix qm_engine default |
| `MD/pada1_test_config.yaml` | Add orca_nproc / orca_mem_mb |
| `MD/tests/test_pipeline_orca.py` | Create — unit tests for changed logic |

---

### Task 1: Create test file and write failing tests

**Files:**
- Create: `MD/tests/__init__.py`
- Create: `MD/tests/test_pipeline_orca.py`

- [ ] **Step 1: Create the tests directory and empty init**

```bash
mkdir -p "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD/tests"
touch "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD/tests/__init__.py"
```

- [ ] **Step 2: Write the test file**

Create `MD/tests/test_pipeline_orca.py`:

```python
"""Unit tests for ORCA/MCPB.py integration in pada1_md_pipeline.py."""
import sys
import textwrap
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

# Add MD dir to path so we can import the pipeline
sys.path.insert(0, str(Path(__file__).parent.parent))
import pada1_md_pipeline as pipeline


# ── DEFAULT_CONFIG tests ──────────────────────────────────────────────────────

def test_default_qm_engine_is_orca():
    assert pipeline.DEFAULT_CONFIG["qm_engine"] == "orca"

def test_default_config_has_orca_exe():
    assert "orca_exe" in pipeline.DEFAULT_CONFIG

def test_default_config_has_orca_nproc():
    assert "orca_nproc" in pipeline.DEFAULT_CONFIG
    assert isinstance(pipeline.DEFAULT_CONFIG["orca_nproc"], int)

def test_default_config_has_orca_mem_mb():
    assert "orca_mem_mb" in pipeline.DEFAULT_CONFIG
    assert isinstance(pipeline.DEFAULT_CONFIG["orca_mem_mb"], int)

def test_default_config_no_gamess_keys():
    assert "gamess_exe" not in pipeline.DEFAULT_CONFIG
    assert "gamess_version" not in pipeline.DEFAULT_CONFIG


# ── stage2: MCPB.py input content ────────────────────────────────────────────

def _make_minimal_merged_pdb(path: Path):
    """Write a minimal PDB with one Fe HETATM so stage2 can find it."""
    path.write_text(textwrap.dedent("""\
        HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00          FE
        END
    """))


def _run_stage2_with_engine(engine: str, tmp_path: Path) -> Path:
    """Call stage2 with a mocked MCPB.py run, return the mcpb .in file path."""
    cfg = pipeline.DEFAULT_CONFIG.copy()
    cfg["workdir"] = str(tmp_path / "work")
    cfg["qm_engine"] = engine
    cfg["orca_nproc"] = 4
    cfg["orca_mem_mb"] = 4000
    cfg["gaussian_nproc"] = 4
    cfg["gaussian_mem"] = "8GB"

    wdir = tmp_path / "work" / "02_mcpb"
    wdir.mkdir(parents=True)

    pdb_id = cfg["pdb_id"].lower()
    merged = wdir / f"{pdb_id}_merged.pdb"
    _make_minimal_merged_pdb(merged)

    # Pre-create the expected QM input files so stage2 doesn't fail on missing files
    ext = "inp" if engine == "orca" else "com"
    for stem in ["small_opt", "small_fc", "large_mk"]:
        f = wdir / f"{pdb_id}_heme_{stem}.{ext}"
        f.write_text("# placeholder\n0 6\nFe 0 0 0\n")

    # Also need prep outputs for stage2
    prep_dir = tmp_path / "work" / "01_prep"
    prep_dir.mkdir(parents=True)
    for name in [f"{cfg['pdb_id']}_protein.pdb", f"{cfg['pdb_id']}_heme.pdb", f"{cfg['pdb_id']}_mg.pdb"]:
        (prep_dir / name).write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n")

    paths = {
        "prot_pdb": str(prep_dir / f"{cfg['pdb_id']}_protein.pdb"),
        "heme_pdb": str(prep_dir / f"{cfg['pdb_id']}_heme.pdb"),
    }

    with patch("pada1_md_pipeline.run"), \
         patch("pada1_md_pipeline.require_binary"):
        pipeline.stage2_mcpb_step1(cfg, paths)

    return wdir / f"{pdb_id}_heme.in"


def test_mcpb_input_contains_quantum_soft_orca(tmp_path):
    mcpb_in = _run_stage2_with_engine("orca", tmp_path)
    content = mcpb_in.read_text()
    assert "quantum_soft orca" in content


def test_mcpb_input_contains_orca_mem_mb(tmp_path):
    mcpb_in = _run_stage2_with_engine("orca", tmp_path)
    content = mcpb_in.read_text()
    assert "orca_mem_mb" in content
    assert "4000" in content


def test_mcpb_input_no_quantum_soft_for_gaussian(tmp_path):
    mcpb_in = _run_stage2_with_engine("gaussian", tmp_path)
    content = mcpb_in.read_text()
    assert "quantum_soft" not in content


# ── stage2: output file extensions ───────────────────────────────────────────

def _run_stage2_get_paths(engine: str, tmp_path: Path) -> dict:
    cfg = pipeline.DEFAULT_CONFIG.copy()
    cfg["workdir"] = str(tmp_path / "work")
    cfg["qm_engine"] = engine
    cfg["orca_nproc"] = 4
    cfg["orca_mem_mb"] = 4000
    cfg["gaussian_nproc"] = 4
    cfg["gaussian_mem"] = "8GB"

    wdir = tmp_path / "work" / "02_mcpb"
    wdir.mkdir(parents=True)

    pdb_id = cfg["pdb_id"].lower()
    merged = wdir / f"{pdb_id}_merged.pdb"
    _make_minimal_merged_pdb(merged)

    ext = "inp" if engine == "orca" else "com"
    for stem in ["small_opt", "small_fc", "large_mk"]:
        (wdir / f"{pdb_id}_heme_{stem}.{ext}").write_text("# placeholder\n0 6\n")

    prep_dir = tmp_path / "work" / "01_prep"
    prep_dir.mkdir(parents=True)
    paths = {
        "prot_pdb": str(prep_dir / f"{cfg['pdb_id']}_protein.pdb"),
        "heme_pdb": str(prep_dir / f"{cfg['pdb_id']}_heme.pdb"),
    }
    for p in paths.values():
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        Path(p).write_text("ATOM      1  CA  ALA A   1\nEND\n")

    with patch("pada1_md_pipeline.run"), \
         patch("pada1_md_pipeline.require_binary"):
        return pipeline.stage2_mcpb_step1(cfg, paths)


def test_stage2_orca_returns_inp_paths(tmp_path):
    paths = _run_stage2_get_paths("orca", tmp_path)
    assert paths["small_opt"].endswith(".inp")
    assert paths["small_fc"].endswith(".inp")
    assert paths["large_mk"].endswith(".inp")


def test_stage2_gaussian_returns_com_paths(tmp_path):
    paths = _run_stage2_get_paths("gaussian", tmp_path)
    assert paths["small_opt"].endswith(".com")
    assert paths["small_fc"].endswith(".com")
    assert paths["large_mk"].endswith(".com")


# ── _recover_paths: file extension by engine ─────────────────────────────────

def test_recover_paths_orca_uses_inp(tmp_path):
    cfg = pipeline.DEFAULT_CONFIG.copy()
    cfg["workdir"] = str(tmp_path)
    cfg["qm_engine"] = "orca"
    pdb_id = cfg["pdb_id"].lower()

    # Create the .inp files so they register as found
    mcpb_dir = tmp_path / "02_mcpb"
    mcpb_dir.mkdir()
    for stem in ["small_opt", "small_fc", "large_mk"]:
        (mcpb_dir / f"{pdb_id}_heme_{stem}.inp").write_text("")

    paths = pipeline._recover_paths(cfg, start_stage=4)
    assert paths["small_opt"].endswith(".inp")
    assert paths["small_fc"].endswith(".inp")
    assert paths["large_mk"].endswith(".inp")


def test_recover_paths_gaussian_uses_com(tmp_path):
    cfg = pipeline.DEFAULT_CONFIG.copy()
    cfg["workdir"] = str(tmp_path)
    cfg["qm_engine"] = "gaussian"
    pdb_id = cfg["pdb_id"].lower()

    mcpb_dir = tmp_path / "02_mcpb"
    mcpb_dir.mkdir()
    for stem in ["small_opt", "small_fc", "large_mk"]:
        (mcpb_dir / f"{pdb_id}_heme_{stem}.com").write_text("")

    paths = pipeline._recover_paths(cfg, start_stage=4)
    assert paths["small_opt"].endswith(".com")
    assert paths["small_fc"].endswith(".com")
    assert paths["large_mk"].endswith(".com")
```

- [ ] **Step 3: Run tests to confirm they all fail**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python -m pytest tests/test_pipeline_orca.py -v 2>&1 | head -60
```

Expected: multiple FAILs — `default_qm_engine_is_orca`, `no_gamess_keys`, `quantum_soft_orca`, `inp_paths`, etc.

- [ ] **Step 4: Commit the failing tests**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline"
git add MD/tests/
git commit -m "test: add failing tests for ORCA/MCPB.py integration"
```

---

### Task 2: Fix DEFAULT_CONFIG

**Files:**
- Modify: `MD/pada1_md_pipeline.py:84-92`

- [ ] **Step 1: Replace the QM engine block in DEFAULT_CONFIG**

In `pada1_md_pipeline.py`, replace lines 84–92:

```python
    # ── QM engine (Gaussian or GAMESS-US) ────────────────────────────────────
    "qm_engine":     "gamess",        # "gamess" (free) or "gaussian"
    "gaussian_exe":  "g16",           # used only if qm_engine = "gaussian"
    "gamess_exe":    "rungms",        # used only if qm_engine = "gamess"
    "gamess_version": "00",           # GAMESS version string passed to rungms
    "gaussian_nproc": 8,
    "gaussian_mem":  "16GB",
    "spin_mult":     6,               # ferric high-spin S=5/2 → multiplicity=6
    "charge":        0,               # overall QM model charge
```

with:

```python
    # ── QM engine ────────────────────────────────────────────────────────────
    "qm_engine":     "orca",          # "orca" (free, academic) or "gaussian"
    "orca_exe":      "/root/orca_6_1_1/orca",
    "orca_nproc":    4,               # parallel cores for ORCA %pal
    "orca_mem_mb":   4000,            # MB per core for ORCA %maxcore
    "gaussian_exe":  "g16",           # used only if qm_engine = "gaussian"
    "gaussian_nproc": 8,
    "gaussian_mem":  "16GB",
    "spin_mult":     6,               # ferric high-spin S=5/2 → multiplicity=6
    "charge":        0,               # overall QM model charge
```

- [ ] **Step 2: Run the DEFAULT_CONFIG tests**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python -m pytest tests/test_pipeline_orca.py -v -k "default_config or default_qm"
```

Expected: `test_default_qm_engine_is_orca` PASS, `test_default_config_has_orca_exe` PASS, `test_default_config_has_orca_nproc` PASS, `test_default_config_has_orca_mem_mb` PASS, `test_default_config_no_gamess_keys` PASS

- [ ] **Step 3: Commit**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline"
git add MD/pada1_md_pipeline.py
git commit -m "fix: update DEFAULT_CONFIG with ORCA keys, remove GAMESS"
```

---

### Task 3: Fix stage2 — ORCA-aware MCPB.py input and file extensions

**Files:**
- Modify: `MD/pada1_md_pipeline.py:357-393`

- [ ] **Step 1: Move qm_engine detection up and rewrite the MCPB input block**

In `pada1_md_pipeline.py`, replace lines 312–326 (the MCPB.py input write block):

```python
    # ── write MCPB.py input ──────────────────────────────────────────────────
    log.info(f"Writing MCPB.py input → {mcpb_in}")
    mcpb_in.write_text(textwrap.dedent(f"""
        original_pdb {merged.name}
        group_name {pdb_id}_heme
        cut_off {cfg['mcpb_cutoff']}
        ion_ids {fe_serial}
        ion_mol2files FE.mol2
        naa_mol2files HEM.mol2
        frcmod_files HEM.frcmod
        ff_choice {cfg['mcpb_ff']}
        gaff_version {cfg['mcpb_gaff']}
        large_opt 1
        conf_search 0
    """).strip() + "\n")
```

with:

```python
    # ── write MCPB.py input ──────────────────────────────────────────────────
    qm_engine = cfg.get("qm_engine", "orca").lower()
    orca_block = (
        f"\nquantum_soft orca\norca_mem_mb {cfg.get('orca_mem_mb', 4000)}\n"
        if qm_engine == "orca" else ""
    )
    log.info(f"Writing MCPB.py input → {mcpb_in}")
    mcpb_in.write_text(textwrap.dedent(f"""
        original_pdb {merged.name}
        group_name {pdb_id}_heme
        cut_off {cfg['mcpb_cutoff']}
        ion_ids {fe_serial}
        ion_mol2files FE.mol2
        naa_mol2files HEM.mol2
        frcmod_files HEM.frcmod
        ff_choice {cfg['mcpb_ff']}
        gaff_version {cfg['mcpb_gaff']}
        large_opt 1
        conf_search 0
    """).strip() + orca_block + "\n")
```

- [ ] **Step 2: Replace the expected-output and patching block (lines 357–393)**

Replace from `qm_engine = cfg.get(...)` at line 357 through the end of the `else:` block at line 393:

```python
    # MCPB.py generates .inp for ORCA, .com for Gaussian
    ext = "inp" if qm_engine == "orca" else "com"
    log.info(f"Running MCPB.py step 1 (QM engine: {qm_engine})...")
    run(f"MCPB.py -i {mcpb_in.name} -s 1", cwd=str(wdir))

    small_opt = wdir / f"{pdb_id}_heme_small_opt.{ext}"
    small_fc  = wdir / f"{pdb_id}_heme_small_fc.{ext}"
    large_mk  = wdir / f"{pdb_id}_heme_large_mk.{ext}"

    for f in [small_opt, small_fc, large_mk]:
        if not f.exists():
            raise FileNotFoundError(f"Expected QM input not generated: {f}")
        log.info(f"  Generated: {f.name}")

    spin   = cfg["spin_mult"]
    charge = cfg["charge"]

    if qm_engine == "gaussian":
        nproc = cfg["gaussian_nproc"]
        mem   = cfg["gaussian_mem"]
        for com_file in [small_opt, small_fc, large_mk]:
            txt = com_file.read_text()
            txt = re.sub(r"%[Nn]proc\w*=\d+", f"%nproc={nproc}", txt)
            txt = re.sub(r"%[Mm]em=\S+",      f"%mem={mem}",     txt)
            txt = re.sub(r"^(\s*0\s+1\s*)$", f"  {charge} {spin}", txt, flags=re.MULTILINE)
            com_file.write_text(txt)
        log.info(f"  Patched Gaussian inputs: charge={charge}, mult={spin}, nproc={nproc}, mem={mem}")
    # ORCA: MCPB.py writes correct ORCA format directly — no patching needed
```

- [ ] **Step 3: Run stage2 tests**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python -m pytest tests/test_pipeline_orca.py -v -k "mcpb_input or stage2"
```

Expected: all `mcpb_input_*` and `stage2_*` tests PASS

- [ ] **Step 4: Commit**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline"
git add MD/pada1_md_pipeline.py
git commit -m "fix: stage2 MCPB.py input uses quantum_soft orca and .inp extensions"
```

---

### Task 4: Fix _recover_paths — engine-aware file extension

**Files:**
- Modify: `MD/pada1_md_pipeline.py:1063-1083`

- [ ] **Step 1: Add engine-aware extension to _recover_paths**

In `_recover_paths`, insert after line `wd = Path(cfg["workdir"])`:

```python
    qm_ext = "inp" if cfg.get("qm_engine", "orca").lower() == "orca" else "com"
```

Then replace lines 1069–1071 in the `candidates` dict:

```python
        "small_opt":  wd / "02_mcpb" / f"{pdb_id}_heme_small_opt.com",
        "small_fc":   wd / "02_mcpb" / f"{pdb_id}_heme_small_fc.com",
        "large_mk":   wd / "02_mcpb" / f"{pdb_id}_heme_large_mk.com",
```

with:

```python
        "small_opt":  wd / "02_mcpb" / f"{pdb_id}_heme_small_opt.{qm_ext}",
        "small_fc":   wd / "02_mcpb" / f"{pdb_id}_heme_small_fc.{qm_ext}",
        "large_mk":   wd / "02_mcpb" / f"{pdb_id}_heme_large_mk.{qm_ext}",
```

- [ ] **Step 2: Run _recover_paths tests**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python -m pytest tests/test_pipeline_orca.py -v -k "recover"
```

Expected: `test_recover_paths_orca_uses_inp` PASS, `test_recover_paths_gaussian_uses_com` PASS

- [ ] **Step 3: Run all tests to check nothing regressed**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python -m pytest tests/test_pipeline_orca.py -v
```

Expected: all tests PASS

- [ ] **Step 4: Commit**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline"
git add MD/pada1_md_pipeline.py
git commit -m "fix: _recover_paths uses .inp for ORCA, .com for Gaussian"
```

---

### Task 5: Update config YAML files

**Files:**
- Modify: `MD/pada1_config.yaml`
- Modify: `MD/pada1_test_config.yaml`

- [ ] **Step 1: Update pada1_config.yaml**

Replace the QM engine section in `pada1_config.yaml`:

```yaml
# ── QM engine ────────────────────────────────────────────────────────────────
qm_engine:      "orca"
orca_exe:       "/root/orca_6_1_1/orca"
orca_nproc:     4
orca_mem_mb:    4000
gaussian_exe:   "g16"
gaussian_nproc: 8
gaussian_mem:   "16GB"
spin_mult:      6         # ferric high-spin (S=5/2): 2S+1 = 6   ← DO NOT CHANGE
charge:         0         # overall QM model charge
```

Remove the existing `gamess_exe` / `gamess_version` lines if present.

- [ ] **Step 2: Update pada1_test_config.yaml**

Add `orca_nproc` and `orca_mem_mb` after the `orca_exe` line:

```yaml
qm_engine:      "orca"
orca_exe:       "/root/orca_6_1_1/orca"
orca_nproc:     2
orca_mem_mb:    2000
```

- [ ] **Step 3: Verify configs load cleanly**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python -c "
import pada1_md_pipeline as p
cfg = p.load_config('pada1_config.yaml')
print('qm_engine:', cfg['qm_engine'])
print('orca_exe:', cfg['orca_exe'])
print('orca_nproc:', cfg['orca_nproc'])
print('orca_mem_mb:', cfg['orca_mem_mb'])
assert cfg['qm_engine'] == 'orca'
assert 'gamess_exe' not in p.DEFAULT_CONFIG
print('OK')
"
```

Expected output:
```
qm_engine: orca
orca_exe: /root/orca_6_1_1/orca
orca_nproc: 4
orca_mem_mb: 4000
OK
```

- [ ] **Step 4: Commit**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline"
git add MD/pada1_config.yaml MD/pada1_test_config.yaml
git commit -m "config: add orca_nproc/orca_mem_mb, set qm_engine=orca"
```

---

### Task 6: Smoke-test stage 1 end-to-end

Stage 1 (PDB download + prep) has no QM dependency — verify it still runs cleanly with the updated pipeline.

- [ ] **Step 1: Check pdb4amber is available**

```bash
which pdb4amber 2>/dev/null || echo "not found — need AmberTools conda env"
```

If not found, skip to step 3 (document limitation).

- [ ] **Step 2: Run stage 1 only**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline/MD"
python pada1_md_pipeline.py --config pada1_test_config.yaml --stage 1 2>&1 | tail -20
```

Expected: stage 1 completes, outputs `01_prep/5OXU_protein.pdb`, `01_prep/5OXU_heme.pdb`.

If `pdb4amber` not installed, expected failure is `EnvironmentError: 'pdb4amber' not found` — that's correct behaviour, not a regression.

- [ ] **Step 3: Verify ORCA binary is callable**

```bash
echo "! HF STO-3G
* xyz 0 1
H 0 0 0
H 0 0 0.74
*" > /tmp/orca_test.inp
/root/orca_6_1_1/orca /tmp/orca_test.inp 2>&1 | grep -E "ORCA TERMINATED|ERROR" | head -5
```

Expected: `ORCA TERMINATED NORMALLY`

- [ ] **Step 4: Final commit**

```bash
cd "/root/Desktop/DTU/iGEM - Epoxy resin degradation/Pipeline"
git add MD/tests/
git commit -m "test: confirm ORCA binary smoke test passes" --allow-empty
```

---

## Self-Review

**Spec coverage:**
- ✅ DEFAULT_CONFIG: add orca_nproc/orca_mem_mb, remove GAMESS keys → Task 2
- ✅ stage2: quantum_soft orca + .inp extensions → Task 3
- ✅ _stage3_orca: already correct once paths carry .inp → no separate task needed (paths dict from stage2 now has .inp)
- ✅ _recover_paths: engine-aware extension → Task 4
- ✅ Config YAMLs updated → Task 5
- ✅ ORCA binary verification → Task 6

**Placeholder scan:** None found. All steps contain exact code or commands.

**Type consistency:** `qm_ext` used in Task 4 matches pattern from Task 3 (`"inp"` / `"com"`). `orca_mem_mb` key used identically in DEFAULT_CONFIG (Task 2), MCPB input (Task 3), and YAML (Task 5).
