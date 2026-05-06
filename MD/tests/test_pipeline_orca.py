"""Unit tests for ORCA/MCPB.py integration in pada1_md_pipeline.py."""
import sys
import textwrap
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

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


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_minimal_merged_pdb(path: Path):
    """Write a minimal PDB with one Fe HETATM so stage2 can find it."""
    path.write_text(
        "HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00          FE\n"
        "END\n"
    )


def _run_stage2(engine: str, tmp_path: Path) -> tuple[dict, Path]:
    """Call stage2 with mocked MCPB.py run. Returns (paths_dict, mcpb_in_path)."""
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
    _make_minimal_merged_pdb(wdir / f"{pdb_id}_merged.pdb")

    # Pre-create the expected QM input files so stage2 doesn't fail on missing
    ext = "inp" if engine == "orca" else "com"
    for stem in ["small_opt", "small_fc", "large_mk"]:
        (wdir / f"{pdb_id}_heme_{stem}.{ext}").write_text("# placeholder\n0 6\nFe 0 0 0\n")

    prep_dir = tmp_path / "work" / "01_prep"
    prep_dir.mkdir(parents=True)
    paths_in = {
        "prot_pdb": str(prep_dir / f"{cfg['pdb_id']}_protein.pdb"),
        "heme_pdb": str(prep_dir / f"{cfg['pdb_id']}_heme.pdb"),
    }
    Path(paths_in["prot_pdb"]).write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n"
    )
    # heme PDB must contain the Fe atom so stage2 can find it after merging
    Path(paths_in["heme_pdb"]).write_text(
        "HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00          FE\nEND\n"
    )

    with patch("pada1_md_pipeline.run"), patch("pada1_md_pipeline.require_binary"):
        result = pipeline.stage2_mcpb_step1(cfg, paths_in)

    mcpb_in = wdir / f"{pdb_id}_heme.in"
    return result, mcpb_in


# ── stage2: MCPB.py input content ────────────────────────────────────────────

def test_mcpb_input_contains_quantum_soft_orca(tmp_path):
    _, mcpb_in = _run_stage2("orca", tmp_path)
    assert "quantum_soft orca" in mcpb_in.read_text()

def test_mcpb_input_contains_orca_mem_mb(tmp_path):
    _, mcpb_in = _run_stage2("orca", tmp_path)
    content = mcpb_in.read_text()
    assert "orca_mem_mb" in content
    assert "4000" in content

def test_mcpb_input_no_quantum_soft_for_gaussian(tmp_path):
    _, mcpb_in = _run_stage2("gaussian", tmp_path)
    assert "quantum_soft" not in mcpb_in.read_text()


# ── stage2: output file extensions ───────────────────────────────────────────

def test_stage2_orca_returns_inp_paths(tmp_path):
    paths, _ = _run_stage2("orca", tmp_path)
    assert paths["small_opt"].endswith(".inp")
    assert paths["small_fc"].endswith(".inp")
    assert paths["large_mk"].endswith(".inp")

def test_stage2_gaussian_returns_com_paths(tmp_path):
    paths, _ = _run_stage2("gaussian", tmp_path)
    assert paths["small_opt"].endswith(".com")
    assert paths["small_fc"].endswith(".com")
    assert paths["large_mk"].endswith(".com")


# ── _recover_paths: file extension by engine ─────────────────────────────────

def test_recover_paths_orca_uses_inp(tmp_path):
    cfg = pipeline.DEFAULT_CONFIG.copy()
    cfg["workdir"] = str(tmp_path)
    cfg["qm_engine"] = "orca"
    pdb_id = cfg["pdb_id"].lower()

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
