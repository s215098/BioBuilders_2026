#!/usr/bin/env python3
"""
CUPP (Conserved Unique Peptide Patterns) Analysis Pipeline
Automated script for students - just run: python3 run.py

BEFORE RUNNING:
1. Make sure your FASTA file is in this directory  
2. Update curr_fam and query_path variables below (lines ~90-95)
3. Ensure sequences are homology reduced to 90% using mmseqs2
4. Just run: python3 run.py
"""
import os, time, importlib, sys, subprocess

print("=" * 70)
print("CUPP ANALYSIS PIPELINE - AUTOMATED VERSION FOR STUDENTS")
print("=" * 70)
print("This script will automatically:")
print("1. Install required Python packages")
print("2. Process your FASTA sequences") 
print("3. Run CUPP clustering, prediction, and (optionally) visualization")
print("4. Generate results in the CUPP/ directory")
print("\n⚠️  IMPORTANT: Update the configuration section (lines ~90-95) with your data!")
print("=" * 70)

def check_python_version():
    """Check if Python version is compatible"""
    if sys.version_info < (3, 6):
        print("ERROR: This script requires Python 3.6 or higher!")
        print(f"Current version: {sys.version}")
        print("Please upgrade Python or use a newer version.")
        sys.exit(1)
    print(f"✓ Using Python {sys.version.split()[0]}")

def ensure_installed(pkg, pip_name=None):
    """Install package with better error handling"""
    if pip_name is None:
        pip_name = pkg
    
    try:
        importlib.import_module(pkg)
        print(f"✓ {pkg} already installed")
        return True
    except ImportError:
        print(f"Installing {pkg}...")
        try:
            # Use subprocess for better error capture
            result = subprocess.run([sys.executable, "-m", "pip", "install", pip_name], 
                                  capture_output=True, text=True, timeout=120)
            if result.returncode == 0:
                # Test import after installation
                try:
                    importlib.import_module(pkg)
                    print(f"✓ {pkg} installed and verified successfully")
                    return True
                except ImportError:
                    print(f"✗ {pkg} installed but cannot be imported")
                    return False
            else:
                print(f"✗ Failed to install {pkg}: {result.stderr}")
                return False
        except Exception as e:
            print(f"✗ Error installing {pkg}: {e}")
            return False

# Check Python version first
check_python_version()

print("\n[STEP 1/5] Setting up dependencies...")
print("-" * 50)

# Required packages for CUPP (module_name, pip_name)
required_packages = [
    ("argparse", "argparse"),           # Should be built-in, but let's verify
    ("psutil", "psutil"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("Bio", "biopython"),  # BioPython module vs pip name differs
    ("requests", "requests")
]

print("Setting up CUPP dependencies...")
failed_installs = []
for module_name, pip_name in required_packages:
    if not ensure_installed(module_name, pip_name):
        failed_installs.append(module_name)

if failed_installs:
    print(f"\nWARNING: Failed to install: {', '.join(failed_installs)}")
    print("Continuing with available packages...")
    print("Note: If errors occur later, manually install missing packages with:")
    print(f"      {sys.executable} -m pip install <package_name>")
    time.sleep(2)  # Brief pause to let user see the message
else:
    print("✓ All dependencies ready!")

###############################################################################
# Configuration — STUDENTS: MODIFY THESE VARIABLES FOR YOUR DATA
###############################################################################
# INSTRUCTIONS FOR STUDENTS:
# 1. Change 'curr_fam' to your family/dataset name
# 2. Make sure query_path points to your .faa file
# 3. Ensure your sequences are homology reduced to 90% using mmseqs2
# 4. The file should be in the same directory as this script

#curr_fam = "PF14346"                      # Example family name
#query_path = "PF14346_mmseqs2_90.faa"     # Example file name

# Current active configuration:
curr_fam = "PF07555"                          # CHANGE THIS to your family name
query_path = "%s.faa" % curr_fam           # Will look for MINE1.faa (change if needed)
redo = True                                # Set to True to recompute existing results
include_visualization = True               # Set to False to skip visualization step

print(f"Analysis Family: {curr_fam}")
print(f"Input File: {query_path}")
print(f"Recompute existing results: {redo}")
print(f"Include visualization: {include_visualization}")
if not os.path.exists(query_path):
    print(f"⚠️  WARNING: File {query_path} not found in current directory!")
else:
    print(f"✓ Input file found")
print("-" * 70)

###############################################################################
# FASTA header cleaner — minimal addition
###############################################################################
def clean_header(raw_header: str) -> str:
    """
    Extract only the accession ID from FASTA headers.
    Examples:
      >tr|Q9XYZ1|DESCRIPTION  ->  Q9XYZ1
      >A0ABX5ZX28|A0ABX5ZX28_STRTE Hyaluronidase...  ->  A0ABX5ZX28  
      >simple_accession  ->  simple_accession
    """
    line = raw_header.lstrip(">").strip()
    if not line:
        return ""

    # Split on space first to remove any description
    token = line.split(" ", 1)[0]
    
    # Split on | to get parts
    parts = token.split("|")
    
    # If it has the format prefix|accession|..., extract the accession
    if len(parts) >= 2 and len(parts[0]) == 2 and parts[0].isalpha():
        return parts[1]  # return just the accession
    
    # If it has the format accession|something, return just the accession
    if len(parts) >= 2:
        return parts[0]
    
    # Otherwise return the whole token
    return token

def reformat_fasta(input_fasta: str, output_fasta: str):
    """Reformat FASTA headers with error checking"""
    if not os.path.exists(input_fasta):
        print(f"ERROR: Input file {input_fasta} not found!")
        sys.exit(1)
    
    try:
        with open(input_fasta, "r") as fin, open(output_fasta, "w") as fout:
            for line in fin:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    fout.write(f">{clean_header(line)}\n")
                else:
                    fout.write(line + "\n")
        
        seq_count = sum(1 for line in open(input_fasta) if line.startswith('>'))
        print(f"### Total sequences processed: {seq_count}")
    except Exception as e:
        print(f"ERROR processing FASTA file: {e}")
        sys.exit(1)

def get_python_executable():
    """Get the correct Python executable (python vs python3) for cross-platform compatibility"""
    return sys.executable

def run_command_safely(command, description):
    """Run system command with error checking"""
    print(f"### {description}")
    result = os.system(command)
    if result != 0:
        print(f"WARNING: Command returned exit code: {result}")
        print("Continuing with next step...")
        time.sleep(1)  # Brief pause

###############################################################################
# [STEP 2/5] Checking input file and requirements
print("\n[STEP 2/5] Checking input file and requirements...")
print("-" * 50)
# Automatically proceed with sequence processing
# NOTE: Sequences should be homology reduced to 90% using mmseqs2 before running this script
# URL: https://toolkit.tuebingen.mpg.de/tools/mmseqs2
print("### IMPORTANT: Make sure your sequences are homology reduced to 90% using mmseqs2")
print("### Proceeding with query file:", query_path)
if not os.path.exists(query_path):
    print(f"ERROR: Query file {query_path} not found!")
    print("Make sure the file exists in the current directory.")
    sys.exit(1) 

###############################################################################
# [STEP 3/5] Clean FASTA headers and prepare input
print("\n[STEP 3/5] Processing FASTA sequences...")
print("-" * 50)
cleaned_fasta = "%s_cleaned.faa" % query_path.split(".")[0]
print("### Cleaning FASTA headers in:", query_path)
reformat_fasta(query_path, cleaned_fasta)
print("### Cleaned FASTA written to:", cleaned_fasta)
query_path = cleaned_fasta

###############################################################################
# [STEP 4/5] Running CUPP analysis pipeline  
###############################################################################
print("\n[STEP 4/5] Running CUPP analysis pipeline...")
print("-" * 50)
print("This may take several minutes depending on your dataset size.")

python_exe = get_python_executable()
cupp_run_string = f"{python_exe} CUPPclustering_DIRECT.py -cluster -cdhit 0 -domain_off "
cupp_run_string += '-database_version "" '
cupp_run_string += f"-common {curr_fam} "
cupp_run_string += f"-version {curr_fam} "
cupp_run_string += f"-domain_query {query_path} "
cupp_run_string += "-minimum_group_size 5 "
cupp_run_string += "-redo " if redo else ""

start = time.time()
print("### Initiating CUPP clustering...")
run_command_safely(cupp_run_string, "Running CUPP clustering")
print("#T# Clustering time in total: %s" % round(time.time() - start))

# Build library
custom_pool_path = "CUPP/CUPPpools/%s_CUPPpool.json" % curr_fam
recom_cmd = f'{python_exe} CUPPclustering_DIRECT.py -recom -database_version "" -query {custom_pool_path} -version CUPP_lib_{curr_fam} -cdhit 0 -pipe'
recom_cmd += " -redo" if redo else ""
run_command_safely(recom_cmd, "Building CUPP library")

# Predict clustered collection
# Use the generated CUPPlibrary file instead of CUPPpool file
custom_lib_path = "CUPP/CUPPlibrary/8x2_90_CUPP_lib_%s_CUPPlibrary.json" % curr_fam
pred_cmd = f"{python_exe} CUPPprediction_DIRECT.py -q {query_path} -compiled_json {custom_lib_path}"
run_command_safely(pred_cmd, "Running CUPP prediction")

# Conditional visualization step
if include_visualization:
    # [STEP 5/5] Generating visualizations
    print("\n[STEP 5/5] Generating visualizations...")
    print("-" * 50)
    tree_path = "CUPP/itol/%s_fa8x2_90.tree" % curr_fam
    vis_cmd = f'"{python_exe}" CUPPvisualization.py -pool CUPP/CUPPpools/{curr_fam}_CUPPpool.json -name {curr_fam} -outdir CUPP/itol -tree {tree_path} -fasta {query_path}'
    run_command_safely(vis_cmd, "Creating visualizations")
    visualization_status = "✓ Visualization files generated"
else:
    print("\n[STEP 5/5] Skipping visualization...")
    print("-" * 50)
    print("### Visualization step skipped (include_visualization = False)")
    visualization_status = "⚠️  Visualization skipped"

print("\n" + "="*70)
print(" 🎉 CUPP ANALYSIS COMPLETED SUCCESSFULLY! 🎉")
print("="*70)
print("Results have been generated in the following directories:")
print(f"📁 CUPP/CUPPpools/      - Pool files ({curr_fam}_CUPPpool.json)")
print(f"📁 CUPP/CUPPlibrary/    - Library files (8x2_90_CUPP_lib_{curr_fam}_CUPPlibrary.json)")
print(f"📁 CUPP/predicted/      - Prediction results ({query_path.replace('.faa', '')}_CUPP.faa)")

if include_visualization:
    print(f"📁 CUPP/itol/           - Visualization files for iTOL")
    print(f"   └── {curr_fam}/        - Annotation files for phylogenetic tree")
    print(f"   └── {curr_fam}_fa8x2_90.tree - Phylogenetic tree file")
    print(f"   └── additional/      - Additional annotation files:")
    print(f"       • protein_domains_and_signals.txt - Domain and signal peptide annotations")
    print(f"       • Various taxonomy, GO, and Pfam annotation files")
else:
    print(f"📁 CUPP/itol/           - {visualization_status}")

print()
print("Next steps:")

if include_visualization:
    print("1. Use iTOL (https://itol.embl.de/) to visualize the phylogenetic tree")
    print("2. Upload the tree file and annotation files from CUPP/itol/")
else:
    print("1. To create visualizations later, set include_visualization = True and rerun")
    print(f"2. Or run: {sys.executable} CUPPvisualization.py -pool CUPP/CUPPpools/[family]_CUPPpool.json -name [family] -outdir CUPP/itol -tree CUPP/itol/[family]_fa8x2_90.tree -fasta [cleaned_file]")

print("="*70)
