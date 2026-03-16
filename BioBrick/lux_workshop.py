import streamlit as st

st.set_page_config(
    page_title="Operation Bioluminescence",
    page_icon="🧬",
    layout="centered",
)

# ── Styling ───────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;600&display=swap');

html, body, [class*="css"] { font-family: 'DM Sans', sans-serif; }

.stApp { background: #0a0e1a; color: #c8d6e8; }

h1, h2, h3 { font-family: 'Space Mono', monospace !important; }

.step-card {
    background: #0c1220;
    border: 1px solid #1a2a40;
    border-radius: 14px;
    padding: 1.4rem 1.6rem;
    margin-bottom: 1.2rem;
}
.step-card.active { border-color: #2a5a9a; }

.step-number {
    font-family: 'Space Mono', monospace;
    font-size: 0.65rem;
    font-weight: 700;
    letter-spacing: 0.18em;
    text-transform: uppercase;
    color: #2a4a6c;
    margin-bottom: 0.35rem;
}
.step-number.active { color: #4a8ad0; }

.step-title {
    font-family: 'Space Mono', monospace;
    font-size: 1rem;
    font-weight: 700;
    color: #c8d6e8;
    margin-bottom: 0.7rem;
}

.step-body {
    font-size: 0.9rem;
    color: #7a9ab8;
    line-height: 1.75;
}
.step-body strong { color: #b0c8e0; }

.step-link {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    background: #111e36;
    border: 1px solid #1e3a62;
    border-radius: 6px;
    padding: 3px 10px;
    font-family: 'Space Mono', monospace;
    font-size: 0.72rem;
    color: #6aaaf0 !important;
    text-decoration: none !important;
    letter-spacing: 0.04em;
    margin: 3px 2px;
    transition: border-color 0.2s;
}
.step-link:hover { border-color: #4a8ad0; color: #aad4ff !important; }

.code-block {
    background: #060a14;
    border: 1px solid #141e2e;
    border-radius: 8px;
    padding: 0.8rem 1rem;
    font-family: 'Space Mono', monospace;
    font-size: 0.75rem;
    color: #5a9a7a;
    margin: 0.7rem 0;
    white-space: pre;
    overflow-x: auto;
    line-height: 1.6;
}

.sequence-box {
    background: #060a14;
    border: 1px solid #1a3a1a;
    border-left: 3px solid #2a7a2a;
    border-radius: 0 8px 8px 0;
    padding: 0.8rem 1rem;
    font-family: 'Space Mono', monospace;
    font-size: 0.72rem;
    color: #4a9a5a;
    margin: 0.7rem 0;
    word-break: break-all;
    line-height: 1.7;
    letter-spacing: 0.03em;
}

.mol-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 0.8rem;
    margin: 0.8rem 0;
}
.mol-card {
    background: #0a1424;
    border: 1px solid #1a2a40;
    border-radius: 10px;
    padding: 0.9rem 0.8rem;
    text-align: center;
}
.mol-name {
    font-family: 'Space Mono', monospace;
    font-size: 0.78rem;
    font-weight: 700;
    color: #c8d6e8;
    margin-bottom: 0.3rem;
}
.mol-formula {
    font-family: 'Space Mono', monospace;
    font-size: 0.68rem;
    color: #3a6a8a;
    margin-bottom: 0.5rem;
}
.mol-links { display: flex; flex-direction: column; gap: 4px; align-items: center; }

.tag {
    display: inline-block;
    background: #1a2a4a;
    color: #6a9ad0;
    border-radius: 4px;
    padding: 1px 7px;
    font-size: 0.7rem;
    font-family: 'Space Mono', monospace;
    margin: 0 2px;
}

.story-banner {
    background: #0d1a2e;
    border: 1px solid #1e3050;
    border-top: 3px solid #2a6aaa;
    border-radius: 12px;
    padding: 1.6rem 1.8rem;
    margin-bottom: 2rem;
    line-height: 1.8;
    font-size: 0.95rem;
    color: #8aa4c0;
}
.story-banner .label {
    font-family: 'Space Mono', monospace;
    font-size: 0.65rem;
    letter-spacing: 0.2em;
    text-transform: uppercase;
    color: #2a6aaa;
    margin-bottom: 0.6rem;
}

.reaction-eq {
    font-family: 'Space Mono', monospace;
    font-size: 0.76rem;
    color: #3a5a7a;
    text-align: center;
    padding: 0.7rem 1rem;
    border: 1px dashed #1a2a3a;
    border-radius: 8px;
    margin: 1rem 0;
    letter-spacing: 0.04em;
}

.petri-wrap { display: flex; justify-content: center; margin: 1.5rem 0; }
.petri-dish {
    width: 200px; height: 200px; border-radius: 50%;
    border: 3px solid #2a3a5c;
    background: radial-gradient(ellipse at 40% 35%, #0d1f38 0%, #060c18 70%);
    display: flex; align-items: center; justify-content: center;
    box-shadow: 0 0 0 8px #0d1624, 0 0 0 10px #1a2a40;
    transition: all 0.8s ease;
    font-size: 3.2rem;
    flex-direction: column; gap: 6px;
}
.petri-dish.glow-green {
    background: radial-gradient(ellipse at 50% 50%, #00ff88 0%, #00c066 30%, #006633 60%, #001a0d 100%);
    border-color: #00ff88;
    box-shadow: 0 0 0 8px #0a1a10, 0 0 0 10px #00ff4422, 0 0 40px #00ff8899, 0 0 80px #00ff4444, inset 0 0 40px #00ff8833;
    animation: pulse-green 1.5s ease-in-out infinite;
}
.petri-dish.glow-red {
    background: radial-gradient(ellipse at 50% 50%, #1a0505 0%, #0d0303 70%);
    border-color: #3a1010;
    box-shadow: 0 0 0 8px #0d0808, 0 0 0 10px #1a050522;
}
@keyframes pulse-green {
    0%,100% { box-shadow: 0 0 0 8px #0a1a10, 0 0 0 10px #00ff4422, 0 0 40px #00ff8899, 0 0 80px #00ff4444, inset 0 0 40px #00ff8833; }
    50%     { box-shadow: 0 0 0 8px #0a1a10, 0 0 0 10px #00ff4444, 0 0 60px #00ff88cc, 0 0 120px #00ff4466, inset 0 0 60px #00ff8855; }
}
.petri-label {
    font-family: 'Space Mono', monospace; font-size: 0.62rem;
    color: #2a4a6c; letter-spacing: 0.12em; text-transform: uppercase; margin-top: 0.3rem;
}

.hint-box {
    background: #0d1526; border-left: 3px solid #1e4a8a;
    border-radius: 0 8px 8px 0; padding: 0.8rem 1rem;
    font-size: 0.84rem; color: #7a9ac0; margin: 1rem 0;
    line-height: 1.6;
}

.explain-box {
    background: #080e1c; border: 1px solid #141e30;
    border-radius: 12px; padding: 1.2rem 1.4rem;
    color: #8aa4c0; line-height: 1.75; font-size: 0.9rem; margin-bottom: 1.2rem;
}

.success-msg {
    font-family: 'Space Mono', monospace; font-size: 0.8rem;
    color: #00ee77; text-align: center; letter-spacing: 0.08em; margin-top: 0.5rem;
}
.fail-msg {
    font-family: 'Space Mono', monospace; font-size: 0.8rem;
    color: #aa3333; text-align: center; letter-spacing: 0.08em; margin-top: 0.5rem;
}

div[data-testid="stSelectbox"] label {
    font-family: 'Space Mono', monospace !important; font-size: 0.78rem !important;
    color: #4a7aaa !important; letter-spacing: 0.1em; text-transform: uppercase;
}
div[data-testid="stSelectbox"] > div > div {
    background-color: #0f1829 !important; border-color: #1e2e48 !important;
    color: #c8d6e8 !important; font-family: 'Space Mono', monospace !important;
}
div[data-testid="stNumberInput"] label {
    font-family: 'Space Mono', monospace !important; font-size: 0.75rem !important;
    color: #3a6a9a !important;
}
hr { border-color: #1a2a3a !important; }

.progress-bar-wrap {
    display: flex; gap: 6px; margin-bottom: 2rem; align-items: center;
}
.progress-step {
    flex: 1; height: 4px; border-radius: 2px; background: #1a2a3a; transition: background 0.4s;
}
.progress-step.done   { background: #2a7a4a; }
.progress-step.active { background: #2a5a9a; }
.progress-label {
    font-family: 'Space Mono', monospace; font-size: 0.62rem;
    color: #2a4a6c; letter-spacing: 0.12em; white-space: nowrap;
}
</style>
""", unsafe_allow_html=True)

# ── Session state ─────────────────────────────────────────────────────────────
if "blast_done" not in st.session_state:
    st.session_state.blast_done = False
if "docking_done" not in st.session_state:
    st.session_state.docking_done = False

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("# 🧬 Operation Bioluminescence")
st.markdown(
    "<div style='font-family:Space Mono,monospace;font-size:0.72rem;"
    "color:#2a4a6c;letter-spacing:0.18em;text-transform:uppercase;"
    "margin-bottom:1.8rem;'>iGEM Biobrick Workshop · Protein ID + Molecular Docking</div>",
    unsafe_allow_html=True,
)

# ── Story Banner ──────────────────────────────────────────────────────────────
st.markdown("""
<div class="story-banner">
  <div class="label">// mission briefing</div>
  A biotech startup has engineered <em>E. coli</em> to produce a bioluminescent reporter —
  but the strain stays completely dark. Their lead scientist left behind only one clue:
  a raw protein sequence and a note reading <strong style="color:#c8d6e8;">"find what it eats."</strong>
  <br><br>
  Your team must identify the mystery enzyme using BLAST, locate its 3D crystal structure,
  then use molecular docking to determine which of three candidate substrates actually
  fits the active site. Get it right, and your organism <strong style="color:#00cc66;">lights up</strong>.
  <br><br>
  <span style="color:#3a5a7a;font-size:0.82rem;font-family:'Space Mono',monospace;">The reaction you are trying to unlock:</span>
  <div class="reaction-eq" style="margin-top:0.5rem;">
    LuxAB &nbsp;+&nbsp; FMNH₂ &nbsp;+&nbsp; O₂ &nbsp;+&nbsp; R‑CHO &nbsp;→&nbsp;
    FMN &nbsp;+&nbsp; R‑COOH &nbsp;+&nbsp; H₂O &nbsp;+&nbsp; <strong>hν (490 nm)</strong>
  </div>
</div>
""", unsafe_allow_html=True)

# ── Progress bar ──────────────────────────────────────────────────────────────
s2 = "active" if st.session_state.blast_done else ""
s3 = "active" if st.session_state.docking_done else ""
st.markdown(f"""
<div class="progress-bar-wrap">
  <span class="progress-label">Start</span>
  <div class="progress-step done"></div>
  <div class="progress-step active"></div>
  <div class="progress-step {s2}"></div>
  <div class="progress-step {s3}"></div>
  <span class="progress-label">Glow ✓</span>
</div>
""", unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — BLAST
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card active">
  <div class="step-number active">// Step 01 — BLAST</div>
  <div class="step-title">Identify the mystery enzyme</div>
  <div class="step-body">
    Copy the sequence below and run a <strong>blastp</strong> search against
    UniProtKB/Swiss-Prot to identify the protein, its organism, and its function.
    <br><br>
    <strong>What to record from your results:</strong><br>
    &nbsp;• Top hit name, organism, and gene name<br>
    &nbsp;• % identity and E-value<br>
    &nbsp;• Annotated function and any catalytic residues<br>
    &nbsp;• PDB accession (you'll need it in Step 2)<br><br>
    <a class="step-link" href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch" target="_blank">↗ Open NCBI blastp</a>
    <a class="step-link" href="https://www.uniprot.org" target="_blank">↗ Open UniProt</a>
    <a class="step-link" href="https://www.rcsb.org" target="_blank">↗ Open RCSB PDB</a>
  </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="sequence-box">>mystery_protein_iGEM2026
MSKDIIAKELEELGIPMSSIGWKAGQNIPYSVHTSYFAHHWGDKLTPEDMKQMLETHGIP
VFAHHSYTNNPEYAAAFAEQFGAQVWLSSRTPEDMREMLKQAGIEAQHAPYLTPEQIAEW
LQSHPDVKIPVIGHVWNNPNPEHAASFAEELGAKVWLSSRTPEDMREMLKQAGIEASDAP
YLTPEQIAEWLASHPDVKIPVMGHVWNNPNPEYAASFAQELGAKVWLSSRTPEDMREMLK
QAGIEAQHAPYLTPEQIAEWLASHPDVKIPVMGHVWNNPNPEYAASFAQELGAKVWLSSR
TPDEMQNMLAQHGIPVFAHHSYTNNPEYAAAFAEQFGAQVMGH</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hint-box">
  💡 In NCBI blastp, set the database to <strong>UniProtKB/Swiss-Prot</strong>
  (not nr) for faster, cleaner results. Leave everything else at default.
  Your top hit should have an E-value close to 0.
</div>
""", unsafe_allow_html=True)

blast_confirmed = st.checkbox("✓  I've identified the protein and noted the PDB structure ID", key="blast_cb")
if blast_confirmed:
    st.session_state.blast_done = True
    st.markdown("""
    <div style="font-family:'Space Mono',monospace;font-size:0.78rem;color:#2a8a4a;
    padding:0.4rem 0 0;letter-spacing:0.05em;line-height:1.6;">
    ✓ You should have found: <strong>LuxA</strong> — bacterial luciferase α subunit,
    <em>Aliivibrio fischeri</em>. The best available crystal structure is
    <strong>PDB: 3FGC</strong> (<em>V. harveyi</em>, ~85% identity to <em>A. fischeri</em> LuxA).
    Use 3FGC for docking in the steps below.
    </div>
    """, unsafe_allow_html=True)

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Prepare receptor
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card">
  <div class="step-number">// Step 02 — Receptor prep</div>
  <div class="step-title">Prepare the LuxA structure for docking</div>
  <div class="step-body">
    AutoDock Vina needs the receptor as a <strong>.pdbqt file</strong>
    (PDB format with partial charges and atom types added).
    <br><br>
    <strong>1. Download the structure</strong><br>
    Fetch <span class="tag">3FGC</span> from RCSB as a PDB file.<br><br>
    <a class="step-link" href="https://www.rcsb.org/structure/3FGC" target="_blank">↗ 3FGC on RCSB PDB</a>
    <a class="step-link" href="https://www.cgl.ucsf.edu/chimera/download.html" target="_blank">↗ Download UCSF Chimera</a>
    <a class="step-link" href="https://pymol.org/2/" target="_blank">↗ Download PyMOL</a>
    <br><br>
    <strong>2. Clean and convert in AutoDockTools (ADT)</strong><br>
    &nbsp;• File → Read Molecule → open <code>3FGC.pdb</code><br>
    &nbsp;• Edit → Delete Water<br>
    &nbsp;• Select chain B (LuxB) → delete (keep chain A only)<br>
    &nbsp;• Edit → Hydrogens → Add → Polar Only<br>
    &nbsp;• Edit → Charges → Compute Gasteiger<br>
    &nbsp;• Grid → Macromolecule → Choose → save as <code>3FGC_receptor.pdbqt</code><br><br>
    <a class="step-link" href="https://autodocksuite.scripps.edu/adt/" target="_blank">↗ Download AutoDockTools</a>
    <a class="step-link" href="https://autodock-vina.readthedocs.io/en/latest/docking_basic.html" target="_blank">↗ Vina prep guide</a>
    <br><br>
    <strong>3. Define the docking search box</strong><br>
    The active site is centred on residue <strong>Cys106</strong>.
    Use these coordinates in your Vina config file:
  </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="code-block">center_x = 12.5
center_y =  8.3
center_z = 22.1
size_x   = 20
size_y   = 20
size_z   = 22
exhaustiveness = 8</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hint-box">
  💡 Open <code>3FGC_receptor.pdbqt</code> in Chimera or PyMOL and verify that
  the search box covers the pocket where FMN is co-crystallised — that is the active site.
</div>
""", unsafe_allow_html=True)

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Prepare ligands
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card">
  <div class="step-number">// Step 03 — Ligand prep</div>
  <div class="step-title">Prepare the three candidate substrates</div>
  <div class="step-body">
    Download each molecule as a 3D SDF file from PubChem, then convert to
    <strong>.pdbqt</strong> using Meeko (command line) or AutoDockTools.
    All three molecules share a 12-carbon chain or an aldehyde group —
    but only one is the real LuxA substrate.
  </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="mol-grid">
  <div class="mol-card">
    <div class="mol-name">Dodecanal</div>
    <div class="mol-formula">C₁₂H₂₄O · CID 8194</div>
    <div class="mol-links">
      <a class="step-link" href="https://pubchem.ncbi.nlm.nih.gov/compound/8194" target="_blank">↗ PubChem page</a>
      <a class="step-link" href="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/8194/SDF?record_type=3d" target="_blank">↓ Download 3D SDF</a>
    </div>
  </div>
  <div class="mol-card">
    <div class="mol-name">Hexanal</div>
    <div class="mol-formula">C₆H₁₂O · CID 6184</div>
    <div class="mol-links">
      <a class="step-link" href="https://pubchem.ncbi.nlm.nih.gov/compound/6184" target="_blank">↗ PubChem page</a>
      <a class="step-link" href="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/6184/SDF?record_type=3d" target="_blank">↓ Download 3D SDF</a>
    </div>
  </div>
  <div class="mol-card">
    <div class="mol-name">Dodecanol</div>
    <div class="mol-formula">C₁₂H₂₆O · CID 8193</div>
    <div class="mol-links">
      <a class="step-link" href="https://pubchem.ncbi.nlm.nih.gov/compound/8193" target="_blank">↗ PubChem page</a>
      <a class="step-link" href="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/8193/SDF?record_type=3d" target="_blank">↓ Download 3D SDF</a>
    </div>
  </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hint-box">
  💡 <strong>Converting SDF → .pdbqt with Meeko (recommended):</strong><br><br>
  <span style="color:#4a9a6a;font-family:'Space Mono',monospace;">pip install meeko</span><br>
  <span style="color:#4a9a6a;font-family:'Space Mono',monospace;">mk_prepare_ligand.py -i dodecanal.sdf -o dodecanal.pdbqt</span><br><br>
  Repeat for hexanal and dodecanol. The direct download links above fetch
  the 3D conformer automatically — make sure you don't use a 2D SDF by mistake.
  <br><br>
  <a class="step-link" href="https://github.com/forlilab/Meeko" target="_blank">↗ Meeko on GitHub</a>
</div>
""", unsafe_allow_html=True)

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Run Vina
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card">
  <div class="step-number">// Step 04 — Docking</div>
  <div class="step-title">Run AutoDock Vina — three times</div>
  <div class="step-body">
    Run one docking job per ligand. Use the <strong>same receptor and same box</strong>
    for all three so your scores are directly comparable.
    <br><br>
    <a class="step-link" href="https://vina.scripps.edu/downloads/" target="_blank">↗ Download AutoDock Vina</a>
    <a class="step-link" href="https://autodock-vina.readthedocs.io/en/latest/" target="_blank">↗ Vina documentation</a>
    <br><br>
    <strong>Create a config file</strong> (<code>config.txt</code>) with your receptor path
    and the box coordinates from Step 2:
  </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="code-block">receptor = 3FGC_receptor.pdbqt
center_x = 12.5
center_y =  8.3
center_z = 22.1
size_x   = 20
size_y   = 20
size_z   = 22
exhaustiveness = 8
num_modes = 5</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hint-box">
  <strong>Run each ligand from the terminal:</strong><br><br>
  <span style="color:#4a9a6a;font-family:'Space Mono',monospace;">
  vina --config config.txt --ligand dodecanal.pdbqt --out dodecanal_out.pdbqt<br>
  vina --config config.txt --ligand hexanal.pdbqt    --out hexanal_out.pdbqt<br>
  vina --config config.txt --ligand dodecanol.pdbqt  --out dodecanol_out.pdbqt
  </span><br><br>
  Each run prints a table of binding modes. Record the <strong>Mode 1</strong>
  affinity (kcal/mol) for each — that is your best score.
  <strong>More negative = better binding.</strong>
</div>
""", unsafe_allow_html=True)

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Answer
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card active">
  <div class="step-number active">// Step 05 — Answer</div>
  <div class="step-title">Enter your results &amp; light it up</div>
  <div class="step-body">
    Enter your best binding affinity (Mode 1, kcal/mol) from each run,
    then select which molecule you believe is the real LuxA substrate.<br><br>
    <em>Hint: the most negative kcal/mol value = strongest binder.</em>
  </div>
</div>
""", unsafe_allow_html=True)

col1, col2, col3 = st.columns(3)
with col1:
    score_a = st.number_input("Dodecanal (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sa")
with col2:
    score_b = st.number_input("Hexanal (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sb")
with col3:
    score_c = st.number_input("Dodecanol (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sc")

st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)

substrate = st.selectbox(
    "Which substrate gives the best binding affinity to LuxA?",
    options=["— select —", "Dodecanal (C12 aldehyde)", "Hexanal (C6 aldehyde)", "Dodecanol (C12 alcohol)"],
)

# ── Petri dish ────────────────────────────────────────────────────────────────
correct  = substrate == "Dodecanal (C12 aldehyde)"
answered = substrate != "— select —"

if answered:
    if correct:
        dish_class, emoji, sub_label = "petri-dish glow-green", "✨", "bioluminescence detected"
        st.session_state.docking_done = True
    else:
        dish_class, emoji, sub_label = "petri-dish glow-red", "🦠", "no light emission"
else:
    dish_class, emoji, sub_label = "petri-dish", "🧫", "awaiting substrate"

st.markdown(f"""
<div class="petri-wrap">
  <div style="text-align:center;">
    <div class="{dish_class}"><span>{emoji}</span></div>
    <div class="petri-label">{sub_label}</div>
  </div>
</div>
""", unsafe_allow_html=True)

# ── Feedback ──────────────────────────────────────────────────────────────────
if answered:
    if correct:
        st.markdown(
            '<div class="success-msg">✓ CORRECT — dodecanal is the natural substrate. Organism is luminescent.</div>',
            unsafe_allow_html=True,
        )
        st.balloons()

        st.markdown("---")
        st.markdown("### Why dodecanal wins")

        st.markdown("""
<div class="explain-box">
  <strong style="color:#c8d6e8;">Chain length matters.</strong>
  The LuxA active site contains a deep hydrophobic channel lined with residues
  Leu42, Val77, and Leu109 that snugly accommodates a long aliphatic tail (~C10–C14).
  Dodecanal (C12) fills this channel completely and positions its aldehyde carbonyl
  directly adjacent to the catalytic <strong style="color:#c8d6e8;">Cys106</strong>.<br><br>
  <strong style="color:#c8d6e8;">Hexanal (C6)</strong> has the right functional group
  but the short chain leaves the hydrophobic pocket only half-occupied —
  resulting in weaker binding and poor catalytic geometry.<br><br>
  <strong style="color:#c8d6e8;">Dodecanol (C12)</strong> may actually dock with a
  reasonable score — the chain length is right — but the hydroxyl (–OH) group
  cannot be attacked by Cys106. The mechanism requires an <strong style="color:#c8d6e8;">
  aldehyde carbon (–CHO)</strong> as the electrophilic target.
</div>
""", unsafe_allow_html=True)

        st.markdown("""
<div class="reaction-eq">
Cys106–SH &nbsp;→&nbsp; nucleophilic attack on R‑CHO &nbsp;→&nbsp;
tetrahedral intermediate &nbsp;→&nbsp; chemiluminescent decay &nbsp;→&nbsp; <strong>hν (490 nm)</strong>
</div>
""", unsafe_allow_html=True)

        st.markdown("""
<div class="explain-box" style="margin-top:1rem;">
  <strong style="color:#c8d6e8;">Discussion questions for your team:</strong><br><br>
  1. Your BLAST hit was <em>A. fischeri</em> but you used the <em>V. harveyi</em>
     structure (3FGC). Why? What does this tell you about structural bioinformatics
     in practice?<br><br>
  2. Dodecanol may have docked with a better score than hexanal even though neither
     is catalytically active. What does this reveal about the limits of docking scores
     as a proxy for biological function?<br><br>
  3. LuxA is a classic iGEM biobrick. If you wanted your strain to glow without
     adding exogenous aldehyde, which other <em>lux</em> genes would you co-express,
     and why?
</div>
""", unsafe_allow_html=True)

    else:
        st.markdown(
            '<div class="fail-msg">✗ No bioluminescence. Revisit your scores — the most negative kcal/mol wins.</div>',
            unsafe_allow_html=True,
        )
        if substrate == "Hexanal (C6 aldehyde)":
            st.markdown("""
<div class="hint-box">
  Hexanal carries the right aldehyde group, but the hydrophobic channel of LuxA
  is sized for a C10–C14 chain — a C6 tail leaves half the pocket empty.
  Compare your hexanal score with your dodecanal score. Which is more negative?
</div>
""", unsafe_allow_html=True)
        elif substrate == "Dodecanol (C12 alcohol)":
            st.markdown("""
<div class="hint-box">
  Dodecanol fits the hydrophobic pocket well and your score may even look reasonable —
  but check the functional group. LuxA's catalytic Cys106 needs an
  <strong>aldehyde carbon (–CHO)</strong> to attack. A hydroxyl (–OH) cannot
  form the required tetrahedral intermediate. Look at your dodecanal score again.
</div>
""", unsafe_allow_html=True)