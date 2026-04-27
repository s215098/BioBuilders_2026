#!/usr/bin/env python3
"""
UPO PaDa-I × NNBT  —  AI-Rational Iterative Pipeline
======================================================
Stack: Boltz-2 (co-folding + affinity) + Claude Code CLI (reasoning)

Substrate: N,N-Bis(2-hydroxypropyl)-p-toluidine (NNBT)

Each round:
  1. Write YAML  →  boltz predict  →  parse affinity + structure
  2. Extract contact residues from predicted CIF
  3. Send full context to Claude → get structured mutation rationale
  4. Apply best mutation to sequence string  →  next round

No GNINA. No FoldX. No PDB manipulation.
Mutations are pure Python string operations on the protein sequence.

Install
-------
    pip install boltz biopython requests pyyaml pandas numpy

Download Boltz-2 weights (first run auto-downloads, or:)
    boltz predict --help

Usage
-----
    python pipeline.py --rounds 4 --workdir results/

    # with pocket constraints (recommended):
    python pipeline.py --rounds 4 \
        --pocket_residues 69,121,199,76,191,316
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
HEME_CHAIN = "C"
BOLTZ_RESULTS_PREFIX = "boltz_results_"

# PaDa-I mature sequence (post signal peptide, from 5OXU / Püllmann et al.)
# Starts with Ala (AGCA Golden Gate overhang convention)
# REPLACE with the verified sequence from your Excel supplementary file
PADA1_SEQUENCE = (
    "EPGLPPGPLENSSAKLVNDEAHPWKPLRPGDIRGPCPGLNTLASHGYLPRNGVATPAQIINAVQEGFNFDNQAAIFATYAA"
    "HLVDGNLITDLLSIGRKTRLTGPDPPPPASVGGLNEHGTFEGDASMTRGDAFFGNNHDFNETLFEQLVDYSNRFGGGKYNL"
    "TVAGELRFKRIQDSIATNPNFSFVDFRFFTAYGETTFPANLFVDGRRDDGQLDMDAARSFFQFSRMPDDFFRAPSPRSGTG"
    "VEVVVQAHPMQPGRNVGKINSYTVDPTSSDFSTPCLMYEKFVNITVKSLYPNPTVQLRKALNTNLDFLFQGVAAGCTQVFP"
    "YGRD"
)

NNBT_SMILES = "CC1=CC=C(C=C1)N(CC(C)O)CC(C)O"  # N,N-Bis(2-hydroxypropyl)-p-toluidine

# Residues that must never be mutated (catalytic essentials from the paper)
PROTECTED_RESIDUES = {36, 189, 196}   # Cys36, Arg189, Glu196

# Known hotspots from Ramirez-Escudero et al. 2018, updated for N-dealkylation of NNBT
# Reaction target: N-dealkylation of N,N-Bis(2-hydroxypropyl)-p-toluidine
# Products: secondary aromatic amine (p-tolyl-NH-CH2CHOHCH3) + lactaldehyde (CH3CHOHCHO)
# Demonstrated by Linde et al. 2022 (doi:10.3389/fctls.2022.883263): PaDa-I gives 12% conversion
# Mechanism: UPO Compound I abstracts the α-C–H from –CH2–CH(OH)–CH3, then C–N bond collapses
HOTSPOT_CONTEXT = """
Target reaction: N-dealkylation of NNBT (N,N-Bis(2-hydroxypropyl)-p-toluidine)
Products: secondary aromatic amine + lactaldehyde (detected by Purpald assay at 530 nm)
Baseline: PaDa-I achieves only 12% conversion — the campaign goal is to improve this.

Catalytic mechanism (critical for geometry reasoning):
  UPO Compound I (Fe(IV)=O, porphyrin π-cation radical) abstracts a hydrogen atom
  from the α-carbon of one hydroxypropyl chain: N–CH2–CH(OH)–CH3.
  The α-carbon radical then undergoes C–N bond homolysis, releasing the alkyl fragment
  as lactaldehyde and leaving the secondary amine on the toluidine nitrogen.
  THEREFORE: the α-carbon (the –CH2– directly bonded to N) must be positioned
  ~3.5–4 Å from the ferryl oxygen (Fe=O), NOT the nitrogen atom itself.
  If the nitrogen faces the iron instead, N-oxidation competes with or replaces
  N-dealkylation — this is a WRONG binding mode for our target reaction.

Key structural residues from the crystal structure analysis (5OXU, PaDa-I):
- F69, F121, F199: Phe triad — the para-tolyl ring of NNBT engages in π–π or
                   T-stacking with these residues. This positions the molecule
                   deep in the channel. Mutations here tune ring depth and tilt,
                   which directly controls whether the α-C or the N faces Fe=O.
                   Prefer aromatic-preserving substitutions (F→Y, F→W) or
                   geometry-adjusting reductions (F→L) over full removal.
- F76, F191: Channel entrance "molecular hinge". The para-tolyl ring (~7 Å wide)
             must pass these residues to enter. L311 (see below) already widened
             this gap to 7.8 Å. Monitor entrance width in predictions — if the
             tolyl ring is blocked at the entrance the binding pose will be wrong.
- E196, R189: Catalytic acid-base pair — NEVER mutate these.
- C36: Axial heme ligand (5th coordination) — NEVER mutate.
- A316: Hotspot in flexible loop G314-G318. A316P (JEd-I) gave 1.5-6x improvement
        in kcat/Km for all substrates. Loop flexibility controls how deep bulky
        aromatics penetrate — high-priority candidate for early rounds.
- V244, F274, P277: Outer entrance hydrophobics — tolyl ring contacts these during
                    entry. Smaller residues here ease passage without disrupting
                    the deep-channel geometry.
- A73, T192: Polar residues flanking the channel. The –OH of the hydroxypropyl
             chain can H-bond here. Productive geometry has one –OH pointing toward
             T192 or A73 while the α-C faces Fe=O. If both –OH groups are buried
             in the hydrophobic core, the binding mode is likely unproductive.
- L311 (was F311): F311L widened the entrance from 4.1 → 7.8 Å — essential for
                   accommodating the tolyl ring. Treat as a baseline reference;
                   if predictions show a narrow entrance, suggest this mutation first.

Substrate geometry (NNBT = N,N-Bis(2-hydroxypropyl)-p-toluidine):
- Tolyl ring (~7 Å wide, planar, rigid) anchors in the Phe triad via π-stacking.
- Tertiary amine N sits between the tolyl ring and the two hydroxypropyl chains.
- Two identical –CH2–CH(OH)–CH3 chains extend from N. Either can be dealkylated.
- The α-C (–CH2–) of the reacting chain must be ~3.5–4 Å from Fe=O.
- The para-methyl group of the tolyl ring points toward the channel entrance.
  A shallow hydrophobic pocket for this methyl (e.g. via A316V) improves anchoring.

The heme channel is 17 Å deep, hydrophobic, cone-shaped (frustum geometry).
"""

SYSTEM_PROMPT = f"""You are an expert computational enzymologist specialising in 
fungal peroxygenases (UPOs), specifically the laboratory-evolved AaeUPO variant 
PaDa-I. You are directing an in silico directed evolution campaign to improve
binding and catalytic activity toward the substrate NNBT
(N,N-Bis(2-hydroxypropyl)-p-toluidine).

{HOTSPOT_CONTEXT}

Your role is to analyse docking/co-folding results from Boltz-2 and propose 
the single most rational point mutation for the next round. You must reason 
from first principles: consider steric effects, hydrophobicity, hydrogen bonding, 
channel geometry, and the known structural biology of this enzyme.

RULES:
- Never propose mutations at positions 36, 189, or 196 (catalytic essentials)
- Prefer conservative substitutions unless there is strong structural justification
- Always explain WHY the mutation should improve N-dealkylation of NNBT specifically
- Be concrete: cite specific interactions, distances, or structural features
- The α-carbon (–CH2– directly bonded to N) must face Fe=O at ~3.5–4 Å for
  productive N-dealkylation. This is the single most important geometry check.
  If the N or the tolyl ring faces the iron instead, the binding mode is wrong
  regardless of affinity score — prioritise fixing geometry over improving affinity
- If the channel entrance is too narrow for the para-tolyl ring, propose widening
  mutations first (F76, F191, building on L311) before optimising deep contacts
- If geometry is correct and affinity is good, suggest mutations that better anchor
  the hydroxypropyl –OH groups via H-bonding (T192, A73) to lock the productive pose
- Never optimise for N-oxidation — the target product is lactaldehyde via C–N cleavage

Respond ONLY with valid JSON matching this schema:
{{
  "reasoning": "detailed structural reasoning (3-5 sentences)",
  "proposed_mutation": {{
    "position": <int>,        
    "from_aa": "<single letter>",
    "to_aa": "<single letter>",
    "mutation_string": "<e.g. A316P>"
  }},
  "expected_effect": "what structural/functional change is expected",
  "confidence": "high|medium|low",
  "alternative_mutations": [
    {{"mutation_string": "X000Y", "rationale": "brief reason"}},
    {{"mutation_string": "X000Y", "rationale": "brief reason"}}
  ],
  "warning": "any concern about this mutation or null"
}}"""


# ─────────────────────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BoltzResult:
    affinity_pred_value: float      # log10(IC50 in µM) — lower = better binder
    affinity_probability: float     # 0-1, probability of being a binder
    confidence_score: float         # pLDDT-like overall confidence
    structure_cif: Path
    contacts: list[str] = field(default_factory=list)


@dataclass
class ClaudeRationale:
    reasoning: str
    proposed_mutation: dict         # {position, from_aa, to_aa, mutation_string}
    expected_effect: str
    confidence: str
    alternatives: list[dict]
    warning: Optional[str]
    raw_response: str


@dataclass
class RoundSummary:
    round_num: int
    sequence: str
    mutations_so_far: list[str]
    boltz_result: BoltzResult
    claude_rationale: ClaudeRationale
    mutation_applied: Optional[str]


# ─────────────────────────────────────────────────────────────────────────────
# Boltz-2 interface
# ─────────────────────────────────────────────────────────────────────────────

class BoltzRunner:

    def __init__(self, workdir: Path, pocket_residues: list[int] = None, template: Path = None):
        self.workdir = workdir
        self.pocket_residues = pocket_residues or [69, 121, 199, 76, 191]
        self.template = template

    def write_yaml(
        self,
        sequence: str,
        round_num,
        label: str = "",
    ) -> Path:
        """Write Boltz-2 YAML input for protein + NNBT with affinity prediction."""

        # Pocket constraints — bias co-folding toward heme channel
        constraints = []
        if self.pocket_residues:
            constraints.append({
                "pocket": {
                    "binder": LIGAND_CHAIN,
                    "contacts": [
                        [PROTEIN_CHAIN, r]
                        for r in self.pocket_residues
                    ]
                }
            })

        config = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": PROTEIN_CHAIN,
                        "sequence": sequence,
                    }
                },
                {
                    "ligand": {
                        "id": LIGAND_CHAIN,
                        "smiles": NNBT_SMILES,
                    }
                },
                {
                    "ligand": {
                        "id": HEME_CHAIN,
                        "ccd": "HEM",
                    }
                }
            ],
            "properties": [
                {"affinity": {"binder": LIGAND_CHAIN}}
            ],
        }

        if constraints:
            config["constraints"] = constraints

        if self.template:
            config["templates"] = [
                {
                    "cif": str(self.template),
                    "chain_id": PROTEIN_CHAIN,
                    "force": True,
                    "threshold": 1.5,
                }
            ]

        slug = f"round{round_num}" + (f"_{label}" if label else "")
        yaml_path = self.workdir / f"{slug}.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
        log.info(f"  Wrote Boltz YAML → {yaml_path}")
        return yaml_path

    def run(self, yaml_path: Path, round_num, label: str = "") -> BoltzResult:
        """Run boltz predict and parse outputs."""
        out_dir = self.workdir / "boltz_out"
        slug = f"round{round_num}" + (f"_{label}" if label else "")

        cmd = [
            "boltz", "predict", str(yaml_path),
            "--out_dir", str(out_dir),
            "--use_msa_server",          # auto-generate MSA via mmseqs2 server
            "--diffusion_samples", "5",  # 5 structural samples, take best
        ]

        log.info(f"Running Boltz-2 ({slug})...")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            log.error(result.stderr[-2000:])
            raise RuntimeError(f"Boltz-2 failed for {slug}")

        return self._parse_output(out_dir, yaml_path.stem)

    def _parse_output(self, out_dir: Path, stem: str) -> BoltzResult:
        """Parse affinity JSON and confidence from Boltz-2 output directory."""
        pred_dir = out_dir / f"{BOLTZ_RESULTS_PREFIX}{stem}" / "predictions" / stem

        all_files = sorted(pred_dir.glob("*"))
        cif_files = [f for f in all_files if f.name.endswith("_model_0.cif")]
        if not cif_files:
            raise FileNotFoundError(f"No CIF output in {pred_dir}")
        best_cif = cif_files[0]

        affinity_files = [f for f in all_files if f.name.startswith("affinity_") and f.name.endswith(".json")]
        affinity_pred = 999.0
        affinity_prob = 0.0
        if affinity_files:
            with open(affinity_files[0]) as f:
                aff = json.load(f)
            # Take best sample (lowest affinity_pred_value)
            samples = aff if isinstance(aff, list) else [aff]
            best = min(samples, key=lambda x: x.get("affinity_pred_value", 999))
            affinity_pred = best.get("affinity_pred_value", 999.0)
            affinity_prob = best.get("affinity_probability_binary", 0.0)
        else:
            log.warning("No affinity output found — check Boltz-2 version")

        conf_files = [f for f in all_files if f.name.startswith("confidence_") and f.name.endswith(".json")]
        confidence = 0.0
        if conf_files:
            with open(conf_files[0]) as f:
                conf = json.load(f)
            confidence = conf.get("confidence_score", 0.0)

        log.info(
            f"  affinity_pred_value = {affinity_pred:.3f} log10(IC50/µM) | "
            f"binder_prob = {affinity_prob:.3f} | "
            f"confidence = {confidence:.3f}"
        )

        # Extract contacts from CIF
        contacts = self._extract_contacts_from_cif(best_cif)

        return BoltzResult(
            affinity_pred_value=affinity_pred,
            affinity_probability=affinity_prob,
            confidence_score=confidence,
            structure_cif=best_cif,
            contacts=contacts,
        )

    @staticmethod
    def _extract_contacts_from_cif(cif_path: Path, cutoff: float = 4.5) -> list[str]:
        """
        Parse mmCIF and find protein residues within cutoff Å of ligand (chain B).
        Returns list like ['F69', 'A316', 'F121']
        """
        try:
            from Bio.PDB import MMCIFParser
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("pred", str(cif_path))
            model = structure[0]

            ligand_atoms = []
            for chain in model:
                if chain.id == LIGAND_CHAIN:
                    for res in chain:
                        ligand_atoms.extend(list(res.get_atoms()))

            if not ligand_atoms:
                log.warning(f"No chain {LIGAND_CHAIN} (ligand) found in CIF")
                return []

            lig_coords = np.array([a.get_vector().get_array() for a in ligand_atoms])

            contacts = set()
            aa_map = {
                "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
                "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
                "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
                "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
            }
            for chain in model:
                if chain.id != PROTEIN_CHAIN:
                    continue
                for res in chain:
                    if res.id[0] != " ":
                        continue
                    resname = res.get_resname()
                    resnum = res.id[1]
                    aa1 = aa_map.get(resname, "X")
                    for atom in res.get_atoms():
                        coord = atom.get_vector().get_array()
                        dists = np.linalg.norm(lig_coords - coord, axis=1)
                        if dists.min() <= cutoff:
                            contacts.add(f"{aa1}{resnum}")
                            break

            contact_list = sorted(contacts, key=lambda x: int(re.sub(r"[^0-9]", "", x)))
            log.info(f"  Contact residues (≤{cutoff}Å): {contact_list}")
            return contact_list

        except Exception as e:
            log.warning(f"Contact extraction failed: {e}")
            return []


# ─────────────────────────────────────────────────────────────────────────────
# Claude reasoning engine
# ─────────────────────────────────────────────────────────────────────────────

class ClaudeReasoner:

    def __init__(self):
        pass

    def reason(
        self,
        round_num: int,
        sequence: str,
        mutations_so_far: list[str],
        boltz_result: BoltzResult,
        history: list,
    ) -> ClaudeRationale:
        """
        Send round context to Claude Code CLI and get structured mutation proposal.
        Full context is included in each call so no API key is required.
        """

        prev_context = ""
        if history:
            prev_context = "\n\nPrevious rounds trajectory:\n"
            prev_aff = None
            for s in history:
                mut = s.claude_rationale.proposed_mutation.get("mutation_string", "none")
                aff = s.boltz_result.affinity_pred_value
                prob = s.boltz_result.affinity_probability
                contacts = ", ".join(s.boltz_result.contacts) if s.boltz_result.contacts else "unknown"
                effect = s.claude_rationale.expected_effect

                if prev_aff is None:
                    delta_str = ""
                else:
                    delta = aff - prev_aff
                    delta_str = f" (Δ{delta:+.4f})"

                if s.mutation_applied:
                    outcome = f"accepted"
                else:
                    conf = s.claude_rationale.confidence
                    outcome = f"rejected, {conf} confidence"

                prev_context += (
                    f"  Round {s.round_num}: {mut} — {effect}\n"
                    f"    affinity {aff:.4f}{delta_str}, binder_prob {prob:.3f}, {outcome}\n"
                    f"    contacts: {contacts}\n"
                )
                prev_aff = aff

        user_message = f"""
ROUND {round_num} RESULTS
{'='*50}

Current sequence mutations applied so far: {mutations_so_far if mutations_so_far else 'none (wildtype PaDa-I)'}

Boltz-2 co-folding results for PaDa-I{'+'.join(mutations_so_far) if mutations_so_far else ''} + NNBT:
  affinity_pred_value : {boltz_result.affinity_pred_value:.4f}  (log10 IC50 µM — lower = tighter binding)
  binder_probability  : {boltz_result.affinity_probability:.4f}  (0→1, probability of being a binder)
  structural confidence: {boltz_result.confidence_score:.4f}

Contact residues within 4.5 Å of NNBT in predicted complex:
  {boltz_result.contacts if boltz_result.contacts else 'extraction failed — reason from known structure'}
{prev_context}

Current sequence (full, for position reference):
{sequence}

Based on these results, analyse:
1. CRITICAL GEOMETRY: Is the α-carbon (–CH2– directly bonded to N) positioned
   ~3.5–4 Å from the heme ferryl oxygen (Fe=O)? If the nitrogen or tolyl ring faces
   the iron instead, state this explicitly — it means the binding mode favours
   N-oxidation over N-dealkylation and must be corrected before anything else.
2. Is the para-tolyl ring engaged in π-stacking with F69/F121/F199, or is the
   channel entrance too narrow and blocking productive insertion?
3. Are the hydroxypropyl –OH groups pointing toward T192/A73 (productive, anchors
   the pose) or buried in the hydrophobic core (unproductive, molecule will slide)?
4. What single point mutation would most improve N-dealkylation efficiency — either
   by correcting α-C geometry, widening the entrance, or locking the productive pose?
5. Are there any red flags (wrong regiochemistry, steric clashes, channel too narrow)?

Propose the next mutation as JSON.
"""

        full_prompt = f"{SYSTEM_PROMPT}\n\n{user_message}"

        log.info(f"  Querying Claude Code CLI for mutation rationale (round {round_num})...")

        result = subprocess.run(
            ["claude", "-p", "--no-session-persistence"],
            input=full_prompt,
            capture_output=True,
            text=True,
            timeout=180,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"Claude CLI failed (exit {result.returncode}): {result.stderr[-1000:]}"
            )

        raw = result.stdout.strip()
        log.info(f"  Claude response received ({len(raw)} chars)")

        return self._parse_rationale(raw)

    @staticmethod
    def _parse_rationale(raw: str) -> ClaudeRationale:
        """Extract JSON from Claude response."""
        # Strip markdown fences if present
        clean = re.sub(r"```json\s*|\s*```", "", raw).strip()
        # Find the JSON object
        match = re.search(r"\{.*\}", clean, re.DOTALL)
        if not match:
            raise ValueError(f"No JSON found in Claude response:\n{raw}")

        data = json.loads(match.group())

        return ClaudeRationale(
            reasoning=data.get("reasoning", ""),
            proposed_mutation=data.get("proposed_mutation", {}),
            expected_effect=data.get("expected_effect", ""),
            confidence=data.get("confidence", "low"),
            alternatives=data.get("alternative_mutations", []),
            warning=data.get("warning"),
            raw_response=raw,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Sequence mutation (pure string operation — no FoldX needed)
# ─────────────────────────────────────────────────────────────────────────────

class SequenceMutator:

    @staticmethod
    def apply(sequence: str, position: int, to_aa: str) -> str:
        """
        Apply point mutation at 1-indexed position.
        Boltz-2 numbers from 1 so we keep 1-based indexing throughout.
        """
        idx = position - 1   # convert to 0-based for Python string
        if idx < 0 or idx >= len(sequence):
            raise ValueError(f"Position {position} out of range (sequence length {len(sequence)})")
        original = sequence[idx]
        mutated = sequence[:idx] + to_aa + sequence[idx+1:]
        log.info(f"  Mutation: {original}{position}{to_aa} applied to sequence")
        return mutated

    @staticmethod
    def validate(sequence: str, position: int, from_aa: str) -> bool:
        """Check that the position actually contains the expected amino acid."""
        idx = position - 1
        actual = sequence[idx]
        if actual != from_aa:
            log.warning(
                f"Position {position}: Claude expected '{from_aa}' "
                f"but sequence has '{actual}' — check numbering offset"
            )
            return False
        return True

    @staticmethod
    def is_protected(position: int) -> bool:
        return position in PROTECTED_RESIDUES


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_mutation_string(mut_str: str) -> tuple[str, int, str]:
    """Parse 'A316P' → ('A', 316, 'P')."""
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', mut_str.strip())
    if not m:
        raise ValueError(f"Cannot parse mutation '{mut_str}' — expected format like 'A316P'")
    return m.group(1), int(m.group(2)), m.group(3)


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

class RationalPipeline:

    def __init__(
        self,
        sequence: str,
        rounds: int,
        workdir: Path,
        pocket_residues: list[int],
        template: Path = None,
    ):
        self.initial_sequence = sequence
        self.rounds = rounds
        self.workdir = workdir
        workdir.mkdir(parents=True, exist_ok=True)

        self.boltz = BoltzRunner(workdir, pocket_residues, template)
        self.claude = ClaudeReasoner()
        self.mutator = SequenceMutator()

        self.history: list[RoundSummary] = []

    def run(self):
        log.info("=" * 60)
        log.info("UPO PaDa-I × NNBT  —  AI-Rational Iterative Pipeline")
        log.info(f"Rounds: {self.rounds}  |  Workdir: {self.workdir}")
        log.info("=" * 60)

        current_sequence = self.initial_sequence
        mutations_applied: list[str] = []
        # Carry forward the mutant result from the previous round to avoid
        # re-running Boltz on the same sequence at the start of the next round.
        carried_result: Optional[BoltzResult] = None

        for rnd in range(1, self.rounds + 1):
            log.info(f"\n{'─'*60}")
            log.info(f"ROUND {rnd}  |  Mutations so far: {mutations_applied or ['none']}")
            log.info(f"{'─'*60}")

            # ── 1. Boltz-2 co-folding + affinity ─────────────────────────────
            if carried_result is not None:
                boltz_result = carried_result
                carried_result = None
                log.info(f"  Reusing Boltz result carried from previous round (no extra run)")
            else:
                yaml_path = self.boltz.write_yaml(current_sequence, rnd)
                boltz_result = self.boltz.run(yaml_path, rnd)

            # ── 2. Claude reasoning ───────────────────────────────────────────
            rationale = self.claude.reason(
                round_num=rnd,
                sequence=current_sequence,
                mutations_so_far=mutations_applied,
                boltz_result=boltz_result,
                history=self.history,
            )

            log.info(f"\n  Claude reasoning:\n  {rationale.reasoning}")
            log.info(f"  Proposed: {rationale.proposed_mutation.get('mutation_string')}")
            log.info(f"  Expected: {rationale.expected_effect}")
            log.info(f"  Confidence: {rationale.confidence}")
            if rationale.warning:
                log.warning(f"  ⚠ Warning: {rationale.warning}")

            # ── 3. Validate and apply mutation ────────────────────────────────
            mutation_applied = None
            prop = rationale.proposed_mutation
            position = prop.get("position")
            from_aa = prop.get("from_aa")
            to_aa = prop.get("to_aa")

            if position and to_aa:
                if self.mutator.is_protected(position):
                    log.warning(
                        f"  Claude proposed mutation at protected position {position} "
                        f"— skipping, will ask for alternative next round"
                    )
                elif self.mutator.validate(current_sequence, position, from_aa):
                    new_sequence = self.mutator.apply(current_sequence, position, to_aa)
                    mut_str = prop.get("mutation_string", f"{from_aa}{position}{to_aa}")

                    # Quick validation run with mutant before committing
                    log.info(f"  Validating mutant {mut_str} with Boltz-2...")
                    yaml_mut = self.boltz.write_yaml(new_sequence, rnd, label=mut_str)
                    boltz_mutant = self.boltz.run(yaml_mut, rnd, label=mut_str)

                    delta = boltz_mutant.affinity_pred_value - boltz_result.affinity_pred_value
                    log.info(
                        f"  Δaffinity_pred_value = {delta:+.4f} "
                        f"({'improved ✓' if delta < 0 else 'worse ✗'})"
                    )

                    if delta < 0 or rationale.confidence == "high":
                        # Accept mutation if affinity improved OR Claude is confident
                        current_sequence = new_sequence
                        mutations_applied.append(mut_str)
                        mutation_applied = mut_str
                        carried_result = boltz_mutant  # reuse in next round — no double run
                        log.info(f"  ✓ Mutation {mut_str} accepted")
                    else:
                        log.info(
                            f"  Mutation {mut_str} not accepted "
                            f"(Δ={delta:+.4f}, confidence={rationale.confidence})"
                        )

            # ── 4. Record round ───────────────────────────────────────────────
            summary = RoundSummary(
                round_num=rnd,
                sequence=current_sequence,
                mutations_so_far=list(mutations_applied),
                boltz_result=boltz_result,
                claude_rationale=rationale,
                mutation_applied=mutation_applied,
            )
            self.history.append(summary)
            self._save_round(summary)

        self._save_final_report()
        log.info("\n✓ Pipeline complete.")
        return self.history

    def run_mutation_scan(
        self,
        mutations: list[str],
        rollback_threshold: float = 0.3,
    ) -> list[dict]:
        """
        Sequential ordered mutation scan with backtracking.

        For each mutation in order:
          - Apply it to the current (best) sequence
          - Run Boltz-2 once
          - If Δaffinity_pred_value ≤ rollback_threshold  →  accept, carry result forward
          - Else  →  skip (revert to previous sequence/result), move to next mutation

        Total Boltz runs = 1 (WT baseline) + len(mutations)  — never more.

        Usage:
            python pipeline.py --scan_mutations A316P,F191L,A73T --rollback_threshold 0.3
        """
        log.info("=" * 60)
        log.info("UPO PaDa-I × NNBT  —  Ordered Mutation Scan")
        log.info(f"Mutations (in order): {mutations}")
        log.info(f"Rollback threshold  : Δ > +{rollback_threshold} log10(IC50) → skip")
        log.info("=" * 60)

        current_sequence = self.initial_sequence
        mutations_applied: list[str] = []
        rows: list[dict] = []

        # ── Baseline (WT or starting sequence) ───────────────────────────────
        yaml_wt = self.boltz.write_yaml(current_sequence, 0, label="WT")
        current_result = self.boltz.run(yaml_wt, 0, label="WT")
        rows.append(self._scan_row(
            step=0, label="WT", mutations=[],
            result=current_result, delta=0.0, accepted=True, reason="baseline",
        ))

        # ── Ordered scan ─────────────────────────────────────────────────────
        for step, mut_str in enumerate(mutations, start=1):
            from_aa, position, to_aa = parse_mutation_string(mut_str)

            if SequenceMutator.is_protected(position):
                log.warning(f"  Step {step}: {mut_str} targets protected position {position} — skipped")
                rows.append(self._scan_row(
                    step=step, label=mut_str, mutations=list(mutations_applied),
                    result=current_result, delta=0.0, accepted=False,
                    reason="protected residue",
                ))
                continue

            if not SequenceMutator.validate(current_sequence, position, from_aa):
                log.warning(f"  Step {step}: {mut_str} — residue mismatch, skipped")
                rows.append(self._scan_row(
                    step=step, label=mut_str, mutations=list(mutations_applied),
                    result=current_result, delta=0.0, accepted=False,
                    reason="residue mismatch",
                ))
                continue

            candidate_sequence = SequenceMutator.apply(current_sequence, position, to_aa)
            cumulative_label = "+".join(mutations_applied + [mut_str])

            log.info(f"\n{'─'*60}")
            log.info(f"STEP {step}  |  Testing {mut_str}  (cumulative: {cumulative_label})")
            log.info(f"{'─'*60}")

            yaml_mut = self.boltz.write_yaml(candidate_sequence, step, label=cumulative_label)
            mut_result = self.boltz.run(yaml_mut, step, label=cumulative_label)

            delta = mut_result.affinity_pred_value - current_result.affinity_pred_value
            log.info(
                f"  Δaffinity_pred_value = {delta:+.4f} "
                f"(threshold = +{rollback_threshold})"
            )

            if delta <= rollback_threshold:
                accepted = True
                reason = f"Δ={delta:+.4f} ≤ threshold"
                current_sequence = candidate_sequence
                current_result = mut_result   # carry forward — no extra Boltz run
                mutations_applied.append(mut_str)
                log.info(f"  ✓ {mut_str} accepted — building on this sequence")
            else:
                accepted = False
                reason = f"Δ={delta:+.4f} > threshold — too detrimental, rolling back"
                log.info(f"  ✗ {mut_str} skipped — {reason}")

            rows.append(self._scan_row(
                step=step, label=cumulative_label, mutations=list(mutations_applied),
                result=mut_result, delta=delta, accepted=accepted, reason=reason,
            ))

            # Save per-step JSON
            step_dir = self.workdir / f"scan_step{step}_{mut_str}"
            step_dir.mkdir(exist_ok=True)
            (step_dir / "step_report.json").write_text(json.dumps(rows[-1], indent=2))

        # ── Summary CSV ───────────────────────────────────────────────────────
        df = pd.DataFrame(rows)
        csv_path = self.workdir / "scan_summary.csv"
        df.to_csv(csv_path, index=False)

        # Final accepted sequence
        seq_path = self.workdir / "scan_final_sequence.txt"
        label = "+".join(mutations_applied) or "WT"
        seq_path.write_text(f">PaDa-I_scan_{label}\n{current_sequence}\n")

        log.info("\n" + "=" * 60)
        log.info("SCAN TRAJECTORY")
        log.info("=" * 60)
        log.info(df[["step","label","accepted","delta_affinity",
                      "affinity_pred_value","binder_probability"]].to_string(index=False))
        log.info(f"\nFinal accepted mutations : {mutations_applied or ['none']}")
        log.info(f"Final sequence           → {seq_path}")
        log.info(f"Summary CSV              → {csv_path}")

        return rows

    @staticmethod
    def _scan_row(
        step: int,
        label: str,
        mutations: list[str],
        result: BoltzResult,
        delta: float,
        accepted: bool,
        reason: str,
    ) -> dict:
        return {
            "step": step,
            "label": label,
            "accepted": accepted,
            "reason": reason,
            "cumulative_mutations": "+".join(mutations) or "WT",
            "delta_affinity": round(delta, 4),
            "affinity_pred_value": round(result.affinity_pred_value, 4),
            "binder_probability": round(result.affinity_probability, 4),
            "confidence_score": round(result.confidence_score, 4),
            "contacts": ", ".join(result.contacts),
        }

    def _save_round(self, s: RoundSummary):
        d = self.workdir / f"round{s.round_num}"
        d.mkdir(exist_ok=True)

        report = {
            "round": s.round_num,
            "mutations_so_far": s.mutations_so_far,
            "mutation_applied_this_round": s.mutation_applied,
            "boltz": {
                "affinity_pred_value": s.boltz_result.affinity_pred_value,
                "affinity_probability": s.boltz_result.affinity_probability,
                "confidence": s.boltz_result.confidence_score,
                "contacts": s.boltz_result.contacts,
            },
            "claude_reasoning": s.claude_rationale.reasoning,
            "claude_proposal": s.claude_rationale.proposed_mutation,
            "expected_effect": s.claude_rationale.expected_effect,
            "claude_confidence": s.claude_rationale.confidence,
            "alternatives": s.claude_rationale.alternatives,
            "warning": s.claude_rationale.warning,
            "full_claude_response": s.claude_rationale.raw_response,
        }
        (d / "round_report.json").write_text(json.dumps(report, indent=2))
        log.info(f"  Saved round report → {d / 'round_report.json'}")

    def _save_final_report(self):
        rows = []
        for s in self.history:
            rows.append({
                "round": s.round_num,
                "mutation_applied": s.mutation_applied or "none",
                "cumulative_mutations": "+".join(s.mutations_so_far) or "WT",
                "affinity_pred_value": round(s.boltz_result.affinity_pred_value, 4),
                "binder_probability": round(s.boltz_result.affinity_probability, 4),
                "confidence": round(s.boltz_result.confidence_score, 4),
                "n_contacts": len(s.boltz_result.contacts),
                "contacts": ", ".join(s.boltz_result.contacts),
                "claude_confidence": s.claude_rationale.confidence,
            })

        df = pd.DataFrame(rows)
        csv_path = self.workdir / "pipeline_summary.csv"
        df.to_csv(csv_path, index=False)

        # Print trajectory
        log.info("\n" + "=" * 60)
        log.info("EVOLUTION TRAJECTORY")
        log.info("=" * 60)
        log.info(df[["round","mutation_applied","cumulative_mutations",
                      "affinity_pred_value","binder_probability"]].to_string(index=False))

        # Save final optimised sequence
        if self.history:
            final_seq = self.history[-1].sequence
            final_muts = self.history[-1].mutations_so_far
            seq_path = self.workdir / "final_sequence.txt"
            seq_path.write_text(
                f">PaDa-I_optimised_{'_'.join(final_muts) or 'WT'}\n{final_seq}\n"
            )
            log.info(f"\nFinal sequence → {seq_path}")
            log.info(f"Mutations: {final_muts or 'none'}")

        log.info(f"Summary CSV → {csv_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="UPO PaDa-I × NNBT AI-rational iterative pipeline"
    )
    p.add_argument("--rounds", type=int, default=4)
    p.add_argument("--workdir", type=Path, default=Path("results"))
    p.add_argument(
        "--pocket_residues", type=str, default="69,121,199,76,191,316",
        help="Comma-separated residue numbers for pocket constraint"
    )
    p.add_argument(
        "--sequence", type=str, default=None,
        help="Override default PaDa-I sequence (FASTA string, no header)"
    )
    p.add_argument(
        "--template", type=Path, default=None,
        help="CIF/PDB file to use as backbone template (e.g. 5OXU.cif); anchors chain A at 1.5 Å"
    )
    p.add_argument(
        "--scan_mutations", type=str, default=None,
        help=(
            "Comma-separated ordered mutations to scan (e.g. 'A316P,F191L,A73T'). "
            "Activates scan mode instead of the iterative AI pipeline. "
            "Each mutation is tested in order; if it worsens affinity beyond "
            "--rollback_threshold it is skipped and the scan continues from the "
            "previous best sequence."
        ),
    )
    p.add_argument(
        "--rollback_threshold", type=float, default=0.3,
        help=(
            "Max allowed Δaffinity_pred_value (log10 IC50) for a mutation to be "
            "accepted in scan mode. A mutation that worsens affinity by more than "
            "this value is skipped. Default: 0.3 (≈ 2× worse IC50)."
        ),
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    pocket = [int(x) for x in args.pocket_residues.split(",")]
    sequence = args.sequence or PADA1_SEQUENCE

    pipeline = RationalPipeline(
        sequence=sequence,
        rounds=args.rounds,
        workdir=args.workdir,
        pocket_residues=pocket,
        template=args.template,
    )

    if args.scan_mutations:
        mutations = [m.strip() for m in args.scan_mutations.split(",") if m.strip()]
        pipeline.run_mutation_scan(mutations, rollback_threshold=args.rollback_threshold)
    else:
        pipeline.run()