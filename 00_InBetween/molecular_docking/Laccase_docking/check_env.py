#!/usr/bin/env python3
"""
check_env.py  –  Verify that all dependencies for the docking pipeline are present.
Run this before laccase_docking_pipeline.py to catch missing tools early.
"""

import shutil
import sys

OK  = "  ✓"
FAIL= "  ✗"
WARN= "  ?"

issues = []

def check(label, ok, fix=""):
    sym = OK if ok else FAIL
    print(f"{sym}  {label}")
    if not ok:
        issues.append((label, fix))

def check_python_pkg(import_name, display_name, fix):
    try:
        __import__(import_name)
        check(display_name, True)
    except ImportError:
        check(display_name, False, fix)

print("\n── Python packages ─────────────────────────────────────────────────")
check_python_pkg("Bio",     "biopython",        "pip install biopython")
check_python_pkg("rdkit",   "rdkit",            "conda install -c conda-forge rdkit  OR  pip install rdkit")
check_python_pkg("numpy",   "numpy",            "pip install numpy")
check_python_pkg("requests","requests",         "pip install requests")

print("\n── PDBQT converters ────────────────────────────────────────────────")
try:
    import meeko  # noqa: F401
    check("meeko (preferred PDBQT converter)", True)
except ImportError:
    check("meeko (preferred PDBQT converter)", False, "pip install meeko")

obabel = shutil.which("obabel")
check("openbabel CLI (fallback converter)", bool(obabel),
      "conda install -c conda-forge openbabel")

print("\n── AutoDock Vina ───────────────────────────────────────────────────")
vina_py = False
try:
    from vina import Vina  # noqa: F401
    vina_py = True
    check("Vina Python bindings", True)
except ImportError:
    check("Vina Python bindings", False,
          "conda install -c conda-forge vina")

vina_cli = shutil.which("vina")
check("Vina CLI binary", bool(vina_cli),
      "conda install -c conda-forge vina  OR  download from https://github.com/ccsb-scripps/AutoDock-Vina/releases")

if not vina_py and not vina_cli:
    issues.append(("Vina (both py and CLI missing)",
                   "conda install -c conda-forge vina"))

print()
if issues:
    print("── Issues to fix ───────────────────────────────────────────────────")
    for name, fix in issues:
        print(f"  {name}")
        print(f"    → {fix}")
    print()
    sys.exit(1)
else:
    print("All dependencies satisfied.  You're ready to run the pipeline.\n")
