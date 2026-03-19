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
    font-size: 0.85rem;
    font-weight: 700;
    letter-spacing: 0.18em;
    text-transform: uppercase;
    color: #2a4a6c;
    margin-bottom: 0.35rem;
}
.step-number.active { color: #4a8ad0; }

.step-title {
    font-family: 'Space Mono', monospace;
    font-size: 1.5rem;
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
    color: #2a6aaa;
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
    "color:#2a6aaa;letter-spacing:0.18em;text-transform:uppercase;"
    "margin-bottom:1.8rem;'>iGEM Biobrick Workshop · Protein ID + Molecular Docking</div>",
    unsafe_allow_html=True,
)

# ── Story Banner ──────────────────────────────────────────────────────────────
st.markdown("""
<div class="story-banner">
  <div class="label">// mission briefing</div>
  A biotech startup has engineered <em>E. coli</em> to produce a bioluminescent reporter —
  but the strain stays completely dark. Their lead scientist is suddenly missing and left behind only one clue for them to fix it:
  a raw protein sequence and a note reading <strong style="color:#c8d6e8;">"find what it eats."</strong>
  <br><br>
  Your team must identify the mystery enzyme using BLAST, locate its 3D crystal structure,
  then use molecular docking to determine which of three candidate substrates actually
  fits the active site. Get it right, and your organism <strong style="color:#00cc66;">lights up</strong>.
  <br><br>
  <span style="color:#2a6aaa;font-size:0.82rem;font-family:'Space Mono',monospace;">The reaction you are trying to unlock:</span>
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
<div class='step-card active'>
  <div class='step-number active'>// Step 01 — Identifying Mystery Enzyme</div>
  <div class='step-title'>Identify the mystery enzyme</div>
  <strong>1) NCBI </strong>
  <div class='step-body'>
    Copy the sequence below and run a <strong>blastp</strong> search against
    UniProtKB/Swiss-Prot to identify the protein, its organism, and its function.
    <a class='step-link' href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch" target="_blank">↗ Open NCBI blastp</a>
  </div>
  <div class='sequence-box'>>mystery_protein_iGEM2026
    MEDAKNIKKGPAPFYPLEDGTAGEQLHKAMKRYALVPGTIAFTDAHIEVNITYAEYFEMSVRLAEAMKRY
    GLNTNHRIVVCSENSLQFFMPVLGALFIGVAVAPANDIYNERELLNSMNISQPTVVFVSKKGLQKILNVQ
    KKLPIIQKIIIMDSKTDYQGFQSMYTFVTSHLPPGFNEYDFVPESFDRDKTIALIMNSSGSTGLPKGVAL
    PHRTACVRFSHARDPIFGNQIIPDTAILSVVPFHHGFGMFTTLGYLICGFRVVLMYRFEEELFLRSLQDY
    KIQSALLVPTLFSFFAKSTLIDKYDLSNLHEIASGGAPLSKEVGEAVAKRFHLPGIRQGYGLTETTSAIL
    ITPEGDDKPGAVGKVVPFFEAKVVDLDTGKTLGVNQRGELCVRGPMIMSGYVNNPEATNALIDKDGWLHS
    GDIAYWDEDEHFFIVDRLKSLIKYKGYQVAPAELESILLQHPNIFDAGVAGLPDDDAGELPAAVVVLEHG
    KTMTEKEIVDYVASQVTTAKKLRGGVVFVDEVPKGLTGKLDARKIREILIKAKKGGKSKL</div>
    <div class='step-body'>
    Note down the accession number of the best hit (100% identity) for the next step.
  </div>
  <br>
            
  <strong>2) Uniprot </strong>
    <div class='step-body'>
    go to Uniprot and search for the protein structure.
    ...
    Navigate to the structure section and from the first three hits pick the one with the highest resolution (lowest Å).
    Note down the PDB ID for the next step.
    <a class='step-link' href="https://www.uniprot.org" target="_blank">↗ Open UniProt</a>
  </div>
  <br>

  <strong>2) PDB </strong>
    <div class='step-body'>
    Open the PDB and search for it here. Have a quick look at the PDB entry and structure.      
    <a class='step-link' href="https://www.rcsb.org" target="_blank">↗ Open RCSB PDB</a>
  </div>
</div>
</div>
    
""", unsafe_allow_html=True)

with st.expander("💡 Hint"):
    st.markdown("""
    <div class='hint-box'>
      Your hint text here.
    </div>
    """, unsafe_allow_html=True)


st.markdown("---")




# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Prepare receptor
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card">
  <div class="step-number active">// Step 02 — Mystery Substrates </div>
  <div class="step-title">Download the substrates</div>
  <div class="step-body">
    Download these three mystery substrates. They have already been converted to mol2 format for docking.
  </div>
</div>
""", unsafe_allow_html=True)

# ── Molecule download cards ───────────────────────────────────────────────────
molecules = [
    {
        "name": "Substrate 1",
        "file": "Homophthalic-acid.mol2",
        "label": "↓ Download substrate 1",
    },
    {
        "name": "Substrate 2",
        "file": "2D1S_substrate.mol2",
        "label": "↓ Download substrate 2",
    },
    {
        "name": "Substrate 3",
        "file": "1H-Indole-3-Proprionic_acid.mol2",
        "label": "↓ Download substrate 3",
    },
]

cols = st.columns(3)
for col, mol in zip(cols, molecules):
    with col:
        st.markdown(f"""
        <div class="mol-card">
          <div class="mol-name">{mol["name"]}</div>
        </div>
        """, unsafe_allow_html=True)
        with open(mol["file"], "rb") as f:
            st.download_button(
                label=mol["label"],
                data=f,
                file_name=mol["file"],
                mime="chemical/x-pdbqt",
                use_container_width=True,
                key=mol["file"],
            )

st.markdown("---")





# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Run Vina
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="step-card">
  <div class="step-number active">// Step 03 — Docking</div>
  <div class="step-title">Run AutoDock Vina for each substrate</div>
  <strong>1) Navigate to Autodock Vina </strong>
  <div class="step-body">
    You should now run one docking job with your enzyme and each substrate to be able to compare the scores.
    Open the Swissdock website from the link below. Then navigate to the Autodock Vina tab. 
    <a class="step-link" href="https://www.swissdock.ch" target="_blank">↗ Swissprot website</a>
  </div>
  <strong>2) Upload ligand and target </strong>
  <div class="step-body">
    Upload one of the three ligands (substrate .mol2 files) and press "Prepare Ligand".
    Then upload the target structure (The protein pdb) and press "Prepare Target".
  </div>
  <strong>3) Define search space  </strong>
  <div class="step-body">
    Define the search space for docking by entering the following parameters:
    <br>
    Search box center: 20 - 17 - 8
    <br>
    Search box size: 25 - 30 - 25
    <br>
    You should now see a box appear on the 3D structure — this is the area where the docking algorithm will search for potential binding poses. 
    From literature, we found that the active site should be in here.
  </div>
  <strong>4) Select parameters and run docking </strong>
  <div class="step-body">
    Finally, set the Sampling exhaustivity to 15.
    Then click "Run Docking" and wait for the results to be generated. This may take some minutes and might take longer if many jobs are in the queue.
    When the docking is done, the results appear on the screen.
  </div>
    
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
    score_a = st.number_input("Substrate 1 (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sa")
with col2:
    score_b = st.number_input("Substrate 2 (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sb")
with col3:
    score_c = st.number_input("Substrate 3 (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sc")

# Check the scores:
if st.button("Check Scores"):
    errors = []

    if not (-8 <= score_a <= -4): # The Homophthalic-acid
        errors.append("Substrate 1")
    if not (-10 >= score_b): # The right substrate should be below 10.
        errors.append("Substrate 2")
    if not (-9 <= score_c < -5): # 1H-Indole-3-Proprionic acid
        errors.append("Substrate 3")

    if errors:
        st.warning("One or more scores don't look right — please double check the Mode 1 affinity values from your AutoDock Vina output. Remember: all values should be negative (kcal/mol)." \
        "If you put the right scores already, ask one from the DTU BioBuilder team to help you.")
    else:
        st.success("Scores look good! Now select the correct substrate from the dropdown below and see if your organism glows.")           



st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)

substrate = st.selectbox(
    "Which substrate gives the best binding affinity to your enzyme?",
    options=["— select —", "Substrate 1", "Substrate 2", "Substrate 3"],
)

# ── Petri dish ────────────────────────────────────────────────────────────────
correct  = substrate == "Substrate 2"
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
            '<div class="success-msg">✓ CORRECT! —  Your chosen substrate was the right substrate, DLSA.'
            ' It is a </div>',
            unsafe_allow_html=True,
        )
        st.balloons()

        st.markdown("---")
        st.markdown("### Why DLSA wins")

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
            


# ── References ──────────────────────────────────────────────────────────────────
st.markdown("### 📚 References")
    
st.markdown("**Enzyme & Substrates**")
st.markdown("- Enzyme structure: https://www.rcsb.org/structure/2D1S")
st.markdown("- Substrate 1 - XXX: ")
st.markdown("- Substrate 2 - DLSA: ")
st.markdown("- Substrate 3 - XXX: ")

st.markdown("**Key papers**")
st.markdown("- https://www.nature.com/articles/nature04542")

st.markdown("---")
st.markdown("<div style='font-size:0.75rem; color: gray;'>iGEM Biobrick Workshop</div>", unsafe_allow_html=True)