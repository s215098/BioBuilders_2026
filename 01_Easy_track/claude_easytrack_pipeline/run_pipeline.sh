#!/usr/bin/env bash
# ==============================================================================
#  EasyTrack Drylab Pipeline — Main Runner
#  BioBuilders iGEM 2026
# ==============================================================================
#
#  Usage:
#    bash run_pipeline.sh <config_file> [--from-step N] [--only-step N]
#
#  Examples:
#    # Run the full pipeline with the example config
#    bash run_pipeline.sh config/example_8rnj_badge.yml
#
#    # Resume from step 4 (e.g. you already have sequences + tree)
#    bash run_pipeline.sh config/example_8rnj_badge.yml --from-step 4
#
#    # Run only step 3 (select representatives)
#    bash run_pipeline.sh config/example_8rnj_badge.yml --only-step 3
#
#    # After manually editing selections.txt, apply selections and continue
#    bash run_pipeline.sh config/example_8rnj_badge.yml --from-step 3 --apply-selections
#
# ==============================================================================

set -euo pipefail

# ── Colours ───────────────────────────────────────────────────────────────────
BOLD='\033[1m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
RED='\033[0;31m'
RESET='\033[0m'

# ── Parse arguments ────────────────────────────────────────────────────────────
CONFIG=""
FROM_STEP=1
ONLY_STEP=""
APPLY_SELECTIONS=""

for arg in "$@"; do
  case $arg in
    --from-step)  shift; FROM_STEP="$1"; shift ;;
    --only-step)  shift; ONLY_STEP="$1";  shift ;;
    --apply-selections) APPLY_SELECTIONS="--apply-selections" ;;
    --from-step=*) FROM_STEP="${arg#*=}" ;;
    --only-step=*) ONLY_STEP="${arg#*=}" ;;
    *)
      if [[ -z "$CONFIG" ]]; then
        CONFIG="$arg"
      fi
      ;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo -e "${RED}[ERROR]${RESET} No config file provided."
  echo ""
  echo "  Usage: bash run_pipeline.sh <config_file>"
  echo "  Example: bash run_pipeline.sh config/example_8rnj_badge.yml"
  exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
  echo -e "${RED}[ERROR]${RESET} Config file not found: $CONFIG"
  exit 1
fi

# ── Locate this script (so relative paths work regardless of cwd) ─────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Header ────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${CYAN}============================================================${RESET}"
echo -e "${BOLD}${CYAN}   EasyTrack Drylab Pipeline — BioBuilders iGEM 2026${RESET}"
echo -e "${BOLD}${CYAN}============================================================${RESET}"
echo -e "   Config: ${CONFIG}"
echo -e "   Start from step: ${FROM_STEP}"
if [[ -n "$ONLY_STEP" ]]; then
  echo -e "   Running only step: ${ONLY_STEP}"
fi
echo ""

# ── Dependency check ──────────────────────────────────────────────────────────
echo -e "${BOLD}Checking Python …${RESET}"

# On macOS + conda, the env's Python is 'python' not 'python3'.
# We prefer the conda env's python over the system python3.
if command -v python &>/dev/null && python -c "import yaml" 2>/dev/null; then
  PYTHON=$(command -v python)
elif command -v python3 &>/dev/null && python3 -c "import yaml" 2>/dev/null; then
  PYTHON=$(command -v python3)
else
  echo -e "${RED}[ERROR]${RESET} Could not find a Python with PyYAML."
  echo "        Make sure your conda environment is active:"
  echo "        conda activate easytrack"
  exit 1
fi

echo -e "  ${GREEN}✓${RESET} Python: $($PYTHON --version) ($PYTHON)"

# Warn if the Python looks like the macOS system Python
if [[ "$PYTHON" == "/usr/bin/python"* ]]; then
  echo -e "${YELLOW}  [WARNING]${RESET} Using system Python — conda env may not be active."
  echo "            Run: conda activate easytrack"
fi

echo ""

# ── Helper: run_step ──────────────────────────────────────────────────────────
run_step() {
  local step_num="$1"
  local label="$2"
  local script="$3"
  shift 3
  local extra_args="$@"

  # Skip if before FROM_STEP
  if [[ -n "$ONLY_STEP" && "$step_num" != "$ONLY_STEP" ]]; then
    return 0
  fi
  if [[ "$step_num" -lt "$FROM_STEP" ]]; then
    echo -e "  ${YELLOW}[SKIP]${RESET} Step ${step_num}: ${label}"
    return 0
  fi

  echo ""
  echo -e "${BOLD}${CYAN}──────────────────────────────────────────────────────────${RESET}"
  echo -e "${BOLD}  Step ${step_num}: ${label}${RESET}"
  echo -e "${BOLD}${CYAN}──────────────────────────────────────────────────────────${RESET}"

  $PYTHON "pipeline/${script}" --config "$CONFIG" $extra_args

  local exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    echo ""
    echo -e "${RED}[ERROR]${RESET} Step ${step_num} failed (exit code ${exit_code})."
    echo "        Fix the issue above and re-run with:"
    echo "        bash run_pipeline.sh $CONFIG --from-step $step_num"
    exit $exit_code
  fi

  echo -e "${GREEN}  ✓ Step ${step_num} complete${RESET}"
}

# ── Pipeline steps ─────────────────────────────────────────────────────────────
run_step 1 "Sequence Retrieval"         "01_fetch_sequences.py"
run_step 2 "MSA + Phylogenetic Tree"    "02_build_phylogeny.py"

# Step 3 has a special manual mode: if config says manual, we pause
SELECTION_MODE=$($PYTHON -c "import yaml; cfg=yaml.safe_load(open('$CONFIG')); print(cfg.get('selection',{}).get('mode','auto'))" 2>/dev/null || echo "auto")

if [[ "$SELECTION_MODE" == "manual" && -z "$APPLY_SELECTIONS" ]]; then
  if [[ -n "$ONLY_STEP" && "$ONLY_STEP" != "3" ]]; then
    : # skip
  elif [[ "3" -ge "$FROM_STEP" ]]; then
    echo ""
    echo -e "${BOLD}${CYAN}──────────────────────────────────────────────────────────${RESET}"
    echo -e "${BOLD}  Step 3: Representative Selection  [MANUAL MODE]${RESET}"
    echo -e "${BOLD}${CYAN}──────────────────────────────────────────────────────────${RESET}"
    $PYTHON pipeline/03_select_representatives.py --config "$CONFIG"
    echo ""
    echo -e "${YELLOW}  ⚠  Pipeline paused for manual selection.${RESET}"
    echo -e "     Review the tree image and edit the selections file."
    echo -e "     Then continue with:"
    echo -e "${BOLD}     bash run_pipeline.sh $CONFIG --from-step 3 --apply-selections${RESET}"
    echo ""
    exit 0
  fi
else
  run_step 3 "Representative Selection"   "03_select_representatives.py" $APPLY_SELECTIONS
fi

run_step 4 "Receptor Preparation"       "04_prepare_receptor.py"
run_step 5 "Ligand Preparation"         "05_fetch_ligand.py"
run_step 6 "Molecular Docking"          "06_docking.py"

# Step 6b: Boltz-2 + Targeted Vina (only runs if boltz.enabled: true in config)
BOLTZ_ENABLED=$($PYTHON -c "import yaml; cfg=yaml.safe_load(open('$CONFIG')); print(str(cfg.get('boltz',{}).get('enabled',False)).lower())" 2>/dev/null || echo "false")

if [[ "$BOLTZ_ENABLED" == "true" ]]; then
  run_step "6b" "Boltz-2 + Targeted Vina"  "06b_boltz_vina.py"
else
  if [[ -z "$ONLY_STEP" || "$ONLY_STEP" == "6b" ]]; then
    echo -e "  ${YELLOW}[SKIP]${RESET} Step 6b: Boltz-2 (disabled — set boltz.enabled: true in config to use)"
  fi
fi

run_step 7 "Parse Results"              "07_parse_results.py"

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${GREEN}============================================================${RESET}"
echo -e "${BOLD}${GREEN}   Pipeline complete!${RESET}"
echo -e "${BOLD}${GREEN}============================================================${RESET}"
echo ""

# Print output location from config
OUTPUT_DIR=$($PYTHON -c "import yaml; cfg=yaml.safe_load(open('$CONFIG')); print(cfg.get('output_dir','results'))" 2>/dev/null || echo "results")

echo -e "  Results are in: ${BOLD}${OUTPUT_DIR}/${RESET}"
echo ""
echo -e "  Key output files:"
echo -e "    sequences/representatives.fasta   — curated sequence library"
echo -e "    phylogeny/tree.png                — phylogenetic tree image"
echo -e "    receptor/*_clean_h.pdbqt          — prepared receptor"
echo -e "    ligand/*.sdf                      — ligand structure"
echo -e "    docking/*/poses.pdbqt             — docking poses"
echo -e "    docking/summary.csv               — all binding affinities"
echo -e "    docking/best_poses.txt            — ranked best poses"
echo ""
echo -e "  Visualise in PyMOL:"
echo -e "    pymol ${OUTPUT_DIR}/receptor/*_clean_h.pdb ${OUTPUT_DIR}/docking/*/poses.pdb"
echo ""
