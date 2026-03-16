import streamlit as st
import time

st.set_page_config(
    page_title="LuxA Substrate Finder",
    page_icon="🔬",
    layout="centered",
)

# ── Styling ──────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;600&display=swap');

html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
}

/* Dark lab background */
.stApp {
    background: #0a0e1a;
    color: #c8d6e8;
}

h1, h2, h3 {
    font-family: 'Space Mono', monospace !important;
}

/* Petri dish container */
.petri-wrap {
    display: flex;
    justify-content: center;
    margin: 1.5rem 0;
}

.petri-dish {
    width: 220px;
    height: 220px;
    border-radius: 50%;
    border: 3px solid #2a3a5c;
    background: radial-gradient(ellipse at 40% 35%, #0d1f38 0%, #060c18 70%);
    display: flex;
    align-items: center;
    justify-content: center;
    position: relative;
    box-shadow: 0 0 0 8px #0d1624, 0 0 0 10px #1a2a40;
    transition: all 0.8s ease;
    font-size: 3.5rem;
    flex-direction: column;
    gap: 6px;
}

.petri-dish.glow-green {
    background: radial-gradient(ellipse at 50% 50%, #00ff88 0%, #00c066 30%, #006633 60%, #001a0d 100%);
    border-color: #00ff88;
    box-shadow:
        0 0 0 8px #0a1a10,
        0 0 0 10px #00ff4422,
        0 0 40px #00ff8899,
        0 0 80px #00ff4444,
        inset 0 0 40px #00ff8833;
    animation: pulse-green 1.5s ease-in-out infinite;
}

.petri-dish.glow-red {
    background: radial-gradient(ellipse at 50% 50%, #1a0505 0%, #0d0303 70%);
    border-color: #3a1010;
    box-shadow: 0 0 0 8px #0d0808, 0 0 0 10px #1a050522;
}

@keyframes pulse-green {
    0%, 100% { box-shadow:
        0 0 0 8px #0a1a10,
        0 0 0 10px #00ff4422,
        0 0 40px #00ff8899,
        0 0 80px #00ff4444,
        inset 0 0 40px #00ff8833;
    }
    50% { box-shadow:
        0 0 0 8px #0a1a10,
        0 0 0 10px #00ff4444,
        0 0 60px #00ff88cc,
        0 0 120px #00ff4466,
        inset 0 0 60px #00ff8855;
    }
}

.petri-label {
    font-family: 'Space Mono', monospace;
    font-size: 0.65rem;
    color: #2a4a6c;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-top: 0.3rem;
}

/* Score card styling */
.score-card {
    background: #0f1829;
    border: 1px solid #1e2e48;
    border-radius: 10px;
    padding: 0.9rem 1.2rem;
    font-family: 'Space Mono', monospace;
    font-size: 0.85rem;
    margin: 0.4rem 0;
    display: flex;
    justify-content: space-between;
    align-items: center;
    transition: border-color 0.3s;
}

.score-best {
    border-color: #00cc66;
    color: #00ee77;
}

.score-weak {
    border-color: #2a1a0a;
    color: #664422;
}

.hint-box {
    background: #0d1526;
    border-left: 3px solid #1e4a8a;
    border-radius: 0 8px 8px 0;
    padding: 0.8rem 1rem;
    font-size: 0.85rem;
    color: #7a9ac0;
    margin: 1rem 0;
    font-family: 'Space Mono', monospace;
}

.story-box {
    background: #080e1c;
    border: 1px solid #141e30;
    border-radius: 12px;
    padding: 1.2rem 1.4rem;
    color: #8aa4c0;
    line-height: 1.7;
    font-size: 0.92rem;
    margin-bottom: 1.5rem;
}

.reaction-eq {
    font-family: 'Space Mono', monospace;
    font-size: 0.78rem;
    color: #3a5a7a;
    text-align: center;
    padding: 0.6rem;
    border: 1px dashed #1a2a3a;
    border-radius: 8px;
    margin: 0.8rem 0;
    letter-spacing: 0.05em;
}

.tag {
    display: inline-block;
    background: #1a2a4a;
    color: #6a9ad0;
    border-radius: 4px;
    padding: 1px 7px;
    font-size: 0.72rem;
    font-family: 'Space Mono', monospace;
    letter-spacing: 0.05em;
    margin-right: 4px;
}

/* Override Streamlit select box */
div[data-testid="stSelectbox"] label {
    font-family: 'Space Mono', monospace !important;
    font-size: 0.8rem !important;
    color: #4a7aaa !important;
    letter-spacing: 0.1em;
    text-transform: uppercase;
}

div[data-testid="stSelectbox"] > div > div {
    background-color: #0f1829 !important;
    border-color: #1e2e48 !important;
    color: #c8d6e8 !important;
    font-family: 'Space Mono', monospace !important;
}

.stButton > button {
    background: #0f1829;
    border: 1px solid #1e4a8a;
    color: #6aacf0;
    font-family: 'Space Mono', monospace;
    font-size: 0.8rem;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    border-radius: 6px;
    padding: 0.5rem 1.5rem;
    transition: all 0.2s;
    width: 100%;
}

.stButton > button:hover {
    border-color: #4a8ad0;
    color: #aad4ff;
    background: #141e30;
}

hr {
    border-color: #1a2a3a !important;
}

.success-msg {
    font-family: 'Space Mono', monospace;
    font-size: 0.82rem;
    color: #00ee77;
    text-align: center;
    letter-spacing: 0.08em;
    margin-top: 0.5rem;
}

.fail-msg {
    font-family: 'Space Mono', monospace;
    font-size: 0.82rem;
    color: #aa3333;
    text-align: center;
    letter-spacing: 0.08em;
    margin-top: 0.5rem;
}
</style>
""", unsafe_allow_html=True)


# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("# 🔬 LuxA Substrate Finder")
st.markdown("<div style='font-family: Space Mono, monospace; font-size: 0.75rem; color: #3a5a7a; letter-spacing: 0.12em; text-transform: uppercase; margin-bottom: 1.5rem;'>iGEM Biobrick Workshop · AutoDock Vina Challenge</div>", unsafe_allow_html=True)

# ── Story ─────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="story-box">
<strong style="color: #c8d6e8;">Mission briefing.</strong> A startup has engineered a bacterium that
<em>should</em> glow — but it doesn't. The chassis expresses your mystery enzyme (which you just
identified via BLAST), yet the culture stays dark.<br><br>
You've run AutoDock Vina to dock three candidate substrates into the LuxA active site
(<span class="tag">PDB: 3FGC</span>). Compare your binding energies and select the molecule
that fits best. Choose correctly — and the organism lights up.
</div>
""", unsafe_allow_html=True)

st.markdown('<div class="reaction-eq">LuxAB + FMNH₂ + O₂ + R-CHO  →  FMN + R-COOH + H₂O + <strong>hν (490 nm)</strong></div>', unsafe_allow_html=True)

st.markdown("---")

# ── Score entry ───────────────────────────────────────────────────────────────
st.markdown("### Enter your Vina results")
st.markdown("<div style='font-size:0.82rem; color:#5a7a9a; margin-bottom:1rem;'>Record the best binding affinity (mode 1) from each docking run below, then select your answer.</div>", unsafe_allow_html=True)

col1, col2, col3 = st.columns(3)
with col1:
    score_a = st.number_input("Dodecanal (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sa")
with col2:
    score_b = st.number_input("Hexanal (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sb")
with col3:
    score_c = st.number_input("Dodecanol (kcal/mol)", value=0.0, step=0.1, format="%.1f", key="sc")

# ── Substrate selection ───────────────────────────────────────────────────────
st.markdown("### Your answer")
substrate = st.selectbox(
    "Which substrate gives the best binding affinity to LuxA?",
    options=["— select —", "Dodecanal (C12 aldehyde)", "Hexanal (C6 aldehyde)", "Dodecanol (C12 alcohol)"],
)

# ── Petri dish display ────────────────────────────────────────────────────────
correct = substrate == "Dodecanal (C12 aldehyde)"
answered = substrate != "— select —"

if answered:
    if correct:
        dish_class = "petri-dish glow-green"
        emoji = "✨"
        sub_label = "bioluminescence detected"
    else:
        dish_class = "petri-dish glow-red"
        emoji = "🦠"
        sub_label = "no light emission"
else:
    dish_class = "petri-dish"
    emoji = "🧫"
    sub_label = "awaiting substrate"

st.markdown(f"""
<div class="petri-wrap">
  <div style="text-align:center;">
    <div class="{dish_class}">
      <span>{emoji}</span>
    </div>
    <div class="petri-label">{sub_label}</div>
  </div>
</div>
""", unsafe_allow_html=True)

# ── Feedback ──────────────────────────────────────────────────────────────────
if answered:
    if correct:
        st.markdown('<div class="success-msg">✓ CORRECT — dodecanal is the natural substrate. Organism is luminescent.</div>', unsafe_allow_html=True)
        st.balloons()

        st.markdown("---")
        st.markdown("### Why dodecanal?")
        st.markdown("""
<div class="story-box">
The LuxA active site contains a deep <strong style="color:#c8d6e8;">hydrophobic channel</strong>
lined with residues Leu42, Val77, Leu109 and others that perfectly accommodate a long-chain
aliphatic tail (~C10–C14). Dodecanal (C12) fills this channel and positions its aldehyde
carbonyl close to the catalytic Cys106.<br><br>
<strong style="color:#c8d6e8;">Hexanal</strong> fails because the C6 chain leaves the hydrophobic
pocket only partially occupied — the binding energy penalty is ~2 kcal/mol.<br><br>
<strong style="color:#c8d6e8;">Dodecanol</strong> has the right length but lacks the aldehyde group.
The reaction mechanism requires nucleophilic attack on a carbonyl carbon — a hydroxyl group
cannot form the tetrahedral intermediate needed for light emission.
</div>
""", unsafe_allow_html=True)

        st.markdown('<div class="reaction-eq">Cys106–SH acts as nucleophile → attacks aldehyde C=O → tetrahedral intermediate → chemiluminescent decay → <strong>hν</strong></div>', unsafe_allow_html=True)

    else:
        st.markdown('<div class="fail-msg">✗ No bioluminescence detected. Revisit your docking scores — check which has the most negative kcal/mol.</div>', unsafe_allow_html=True)

        if substrate == "Hexanal (C6 aldehyde)":
            st.markdown("""
<div class="hint-box">
Hint: hexanal has the right functional group (aldehyde) but the hydrophobic pocket of
LuxA is too large for a C6 chain. Compare your score for hexanal vs dodecanal —
which is more negative?
</div>
""", unsafe_allow_html=True)
        elif substrate == "Dodecanol (C12 alcohol)":
            st.markdown("""
<div class="hint-box">
Hint: dodecanol has the right chain length and may even dock reasonably well, but
the -OH group cannot participate in the catalytic mechanism. LuxA requires an
aldehyde (-CHO) for the reaction to proceed. Look at your score for dodecanal.
</div>
""", unsafe_allow_html=True)

# ── Footer info ───────────────────────────────────────────────────────────────
st.markdown("---")
with st.expander("📋 Docking reference info"):
    st.markdown("""
**Receptor:** LuxAB heterodimer from *Vibrio harveyi* — `PDB: 3FGC`
*(closest solved structure to A. fischeri LuxA — ~85% identity)*

**Ligands (download as SDF from PubChem, convert to .pdbqt with AutoDockTools or Meeko):**
| Substrate | PubChem CID | Formula |
|---|---|---|
| Dodecanal | 8194 | C₁₂H₂₄O |
| Hexanal | 6184 | C₆H₁₂O |
| Dodecanol | 8193 | C₁₂H₂₆O |

**Suggested docking box (center on Cys106):**
```
center_x = 12.5
center_y = 8.3
center_z = 22.1
size_x = 20
size_y = 20
size_z = 22
exhaustiveness = 8
```

**Preparing the receptor:**
1. Download `3FGC.pdb` from RCSB
2. Remove chain B (LuxB), waters, and heteroatoms except for FMN if desired
3. Add hydrogens and assign charges in AutoDockTools → save as `3FGC_luxA.pdbqt`
""")

with st.expander("🧬 BLAST exercise hints"):
    st.markdown("""
**The mystery sequence** should be the LuxA subunit from *Aliivibrio fischeri* ES114
(UniProt: Q9KLD5), with 5–8 conservative mutations introduced at non-catalytic surface
residues to make the identity ~92% rather than 100%.

**What students should find:**
- Top hit: LuxA from *A. fischeri* or *V. harveyi* (~85–99% identity)
- Gene name: `luxA`
- Function: bacterial luciferase alpha subunit, flavin monooxygenase
- Key residue: **Cys106** (catalytic), **Arg107** (phosphate binding)
- Next step: retrieve crystal structure → PDB 3FGC

**Discussion question:** Why do we use the *V. harveyi* structure (3FGC) even
though the sequence matched *A. fischeri*? → No *A. fischeri* LuxA crystal structure
is solved at high resolution. This is a real-world bioinformatics decision.
""")