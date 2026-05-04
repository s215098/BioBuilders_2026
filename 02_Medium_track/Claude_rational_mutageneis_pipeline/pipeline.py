#!/usr/bin/env python3
"""
Enzyme × Substrate  —  AI-Rational Iterative Mutation Pipeline
===============================================================
Stack: Boltz-2 (co-folding + affinity) + Claude Code CLI (optional reasoning)

Each round (scan mode):
  1. Write YAML  →  boltz predict  →  parse affinity + contacts from CIF
  2. Accept or roll back mutation based on Δaffinity threshold
  3. If Claude CLI is available: add structural rationale to step report

No GNINA. No FoldX. No PDB manipulation.
Mutations are pure Python string operations on the protein sequence.

Install
-------
    pip install boltz biopython requests pyyaml pandas numpy

Usage
-----
    python pipeline.py --config configs/pada1_nnbt.yaml \\
        --scan_mutations A316P,F191L,A73T \\
        --workdir results/

See configs/pada1_nnbt.yaml for a fully documented example config.
"""

import argparse
import json
import logging
import re
import shutil
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
# Chain IDs (Boltz-2 convention)
# ─────────────────────────────────────────────────────────────────────────────

PROTEIN_CHAIN = "A"
LIGAND_CHAIN  = "B"
BOLTZ_RESULTS_PREFIX = "boltz_results_"


# ─────────────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CofactorConfig:
    chain_id: str
    ccd: Optional[str] = None    # CCD code e.g. "HEM"
    smiles: Optional[str] = None # SMILES string (used if ccd is absent)


@dataclass
class PipelineConfig:
    enzyme_name: str
    sequence: str
    protected_residues: set
    substrate_name: str
    substrate_smiles: str
    pocket_residues: list
    cofactor: Optional[CofactorConfig]
    context: Optional[str]       # free-text reaction/residue context fed to Claude


def load_config(path: Path) -> PipelineConfig:
    with open(path) as f:
        data = yaml.safe_load(f)

    cofactor = None
    if "cofactor" in data and data["cofactor"]:
        c = data["cofactor"]
        cofactor = CofactorConfig(
            chain_id=c["chain_id"],
            ccd=c.get("ccd"),
            smiles=c.get("smiles"),
        )

    return PipelineConfig(
        enzyme_name=data["enzyme"]["name"],
        sequence=data["enzyme"]["sequence"],
        protected_residues=set(data["enzyme"].get("protected_residues", [])),
        substrate_name=data["substrate"]["name"],
        substrate_smiles=data["substrate"]["smiles"],
        pocket_residues=data.get("pocket_residues", []),
        cofactor=cofactor,
        context=data.get("context"),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BoltzResult:
    affinity_pred_value: float   # log10(IC50 in µM) — lower = better binder
    affinity_probability: float  # 0-1, probability of being a binder
    confidence_score: float      # pLDDT-like overall confidence
    structure_cif: Path
    contacts: list = field(default_factory=list)


@dataclass
class ClaudeRationale:
    reasoning: str
    proposed_mutation: dict      # {position, from_aa, to_aa, mutation_string}
    expected_effect: str
    confidence: str
    alternatives: list
    warning: Optional[str]
    raw_response: str


@dataclass
class RoundSummary:
    round_num: int
    sequence: str
    mutations_so_far: list
    boltz_result: BoltzResult
    claude_rationale: Optional[ClaudeRationale]
    mutation_applied: Optional[str]


# ─────────────────────────────────────────────────────────────────────────────
# Boltz-2 interface
# ─────────────────────────────────────────────────────────────────────────────

class BoltzRunner:

    def __init__(self, config: PipelineConfig, workdir: Path, template: Path = None):
        self.config = config
        self.workdir = workdir
        self.template = template

    def write_yaml(self, sequence: str, round_num, label: str = "") -> Path:
        """Write Boltz-2 YAML input for protein + substrate (+ optional cofactor)."""

        constraints = []
        if self.config.pocket_residues:
            constraints.append({
                "pocket": {
                    "binder": LIGAND_CHAIN,
                    "contacts": [
                        [PROTEIN_CHAIN, r]
                        for r in self.config.pocket_residues
                    ]
                }
            })

        sequences = [
            {
                "protein": {
                    "id": PROTEIN_CHAIN,
                    "sequence": sequence,
                }
            },
            {
                "ligand": {
                    "id": LIGAND_CHAIN,
                    "smiles": self.config.substrate_smiles,
                }
            },
        ]

        if self.config.cofactor:
            cof = self.config.cofactor
            entry: dict = {"id": cof.chain_id}
            if cof.ccd:
                entry["ccd"] = cof.ccd
            elif cof.smiles:
                entry["smiles"] = cof.smiles
            sequences.append({"ligand": entry})

        boltz_config = {
            "version": 1,
            "sequences": sequences,
            "properties": [
                {"affinity": {"binder": LIGAND_CHAIN}}
            ],
        }

        if constraints:
            boltz_config["constraints"] = constraints

        if self.template:
            boltz_config["templates"] = [
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
            yaml.dump(boltz_config, f, default_flow_style=False)
        log.info(f"  Wrote Boltz YAML → {yaml_path}")
        return yaml_path

    def run(self, yaml_path: Path, round_num, label: str = "") -> BoltzResult:
        out_dir = self.workdir / "boltz_out"
        slug = f"round{round_num}" + (f"_{label}" if label else "")

        cmd = [
            "boltz", "predict", str(yaml_path),
            "--out_dir", str(out_dir),
            "--use_msa_server",
            "--diffusion_samples", "5",
        ]

        log.info(f"Running Boltz-2 ({slug})...")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            log.error(result.stderr[-2000:])
            raise RuntimeError(f"Boltz-2 failed for {slug}")

        return self._parse_output(out_dir, yaml_path.stem)

    def _parse_output(self, out_dir: Path, stem: str) -> BoltzResult:
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

        contacts = self._extract_contacts_from_cif(best_cif)

        return BoltzResult(
            affinity_pred_value=affinity_pred,
            affinity_probability=affinity_prob,
            confidence_score=confidence,
            structure_cif=best_cif,
            contacts=contacts,
        )

    @staticmethod
    def _extract_contacts_from_cif(cif_path: Path, cutoff: float = 4.5) -> list:
        """Find protein residues within cutoff Å of ligand (chain B)."""
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
# Claude reasoning engine (optional)
# ─────────────────────────────────────────────────────────────────────────────

class ClaudeReasoner:

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._available: Optional[bool] = None

    def is_available(self) -> bool:
        if self._available is None:
            self._available = shutil.which("claude") is not None
        return self._available

    def _build_system_prompt(self) -> str:
        return (
            f"You are an expert computational enzymologist helping engineer "
            f"{self.config.enzyme_name} to improve activity toward "
            f"{self.config.substrate_name}.\n\n"
            f"{self.config.context}\n\n"
            f"Your role is to analyse Boltz-2 co-folding results and propose "
            f"the single most rational point mutation for the next round. "
            f"Reason from first principles: steric effects, hydrophobicity, "
            f"hydrogen bonding, channel geometry, and known structural biology.\n\n"
            f"RULES:\n"
            f"- Never propose mutations at protected positions: "
            f"{sorted(self.config.protected_residues)}\n"
            f"- Prefer conservative substitutions unless there is strong justification\n"
            f"- Always explain WHY the mutation improves activity toward "
            f"{self.config.substrate_name} specifically\n"
            f"- Be concrete: cite specific interactions, distances, or structural features\n\n"
            f"Respond ONLY with valid JSON matching this schema:\n"
            '{{\n'
            '  "reasoning": "detailed structural reasoning (3-5 sentences)",\n'
            '  "proposed_mutation": {{\n'
            '    "position": <int>,\n'
            '    "from_aa": "<single letter>",\n'
            '    "to_aa": "<single letter>",\n'
            '    "mutation_string": "<e.g. A316P>"\n'
            '  }},\n'
            '  "expected_effect": "what structural/functional change is expected",\n'
            '  "confidence": "high|medium|low",\n'
            '  "alternative_mutations": [\n'
            '    {{"mutation_string": "X000Y", "rationale": "brief reason"}}\n'
            '  ],\n'
            '  "warning": "any concern about this mutation or null"\n'
            '}}'
        )

    def reason(
        self,
        round_num: int,
        sequence: str,
        mutations_so_far: list,
        boltz_result: BoltzResult,
        history: list,
    ) -> Optional[ClaudeRationale]:
        """
        Ask Claude for a mutation proposal. Returns None if Claude is not
        available or if the call fails — callers must handle None.
        """
        if not self.is_available():
            log.info("  Claude CLI not found — skipping rationale")
            return None

        if not self.config.context:
            log.warning(
                "  Claude is available but 'context' is not set in config — "
                "skipping rationale. Add a 'context' field to your config YAML."
            )
            return None

        prev_context = ""
        if history:
            prev_context = "\n\nPrevious rounds trajectory:\n"
            prev_aff = None
            for s in history:
                if s.claude_rationale:
                    mut = s.claude_rationale.proposed_mutation.get("mutation_string", "none")
                    effect = s.claude_rationale.expected_effect
                    outcome_conf = s.claude_rationale.confidence
                else:
                    mut = "none"
                    effect = ""
                    outcome_conf = "n/a"
                aff = s.boltz_result.affinity_pred_value
                prob = s.boltz_result.affinity_probability
                contacts = ", ".join(s.boltz_result.contacts) if s.boltz_result.contacts else "unknown"
                delta_str = f" (Δ{aff - prev_aff:+.4f})" if prev_aff is not None else ""
                outcome = "accepted" if s.mutation_applied else f"rejected, {outcome_conf} confidence"
                prev_context += (
                    f"  Round {s.round_num}: {mut} — {effect}\n"
                    f"    affinity {aff:.4f}{delta_str}, binder_prob {prob:.3f}, {outcome}\n"
                    f"    contacts: {contacts}\n"
                )
                prev_aff = aff

        enzyme_label = self.config.enzyme_name
        if mutations_so_far:
            enzyme_label += "+" + "+".join(mutations_so_far)

        user_message = (
            f"ROUND {round_num} RESULTS\n"
            f"{'='*50}\n\n"
            f"Current mutations applied so far: {mutations_so_far or 'none (wildtype)'}\n\n"
            f"Boltz-2 co-folding results for {enzyme_label} + {self.config.substrate_name}:\n"
            f"  affinity_pred_value : {boltz_result.affinity_pred_value:.4f}  (log10 IC50 µM — lower = tighter)\n"
            f"  binder_probability  : {boltz_result.affinity_probability:.4f}  (0→1)\n"
            f"  structural confidence: {boltz_result.confidence_score:.4f}\n\n"
            f"Contact residues within 4.5 Å of {self.config.substrate_name}:\n"
            f"  {boltz_result.contacts or 'extraction failed'}\n"
            f"{prev_context}\n\n"
            f"Current sequence:\n{sequence}\n\n"
            f"Propose the next single point mutation as JSON."
        )

        full_prompt = self._build_system_prompt() + "\n\n" + user_message

        log.info(f"  Querying Claude for mutation rationale (round {round_num})...")
        try:
            result = subprocess.run(
                ["claude", "-p", "--no-session-persistence"],
                input=full_prompt,
                capture_output=True,
                text=True,
                timeout=180,
            )
            if result.returncode != 0:
                log.warning(f"  Claude exited {result.returncode}: {result.stderr[-500:]}")
                return None
        except Exception as e:
            log.warning(f"  Claude call failed: {e}")
            return None

        raw = result.stdout.strip()
        log.info(f"  Claude response received ({len(raw)} chars)")
        return self._parse_rationale(raw)

    @staticmethod
    def _parse_rationale(raw: str) -> Optional[ClaudeRationale]:
        clean = re.sub(r"```json\s*|\s*```", "", raw).strip()
        match = re.search(r"\{.*\}", clean, re.DOTALL)
        if not match:
            log.warning(f"No JSON found in Claude response")
            return None
        try:
            data = json.loads(match.group())
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse Claude JSON: {e}")
            return None

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
# Sequence mutation (pure string operation)
# ─────────────────────────────────────────────────────────────────────────────

class SequenceMutator:

    @staticmethod
    def apply(sequence: str, position: int, to_aa: str) -> str:
        """Apply point mutation at 1-indexed position."""
        idx = position - 1
        if idx < 0 or idx >= len(sequence):
            raise ValueError(f"Position {position} out of range (sequence length {len(sequence)})")
        original = sequence[idx]
        mutated = sequence[:idx] + to_aa + sequence[idx+1:]
        log.info(f"  Mutation: {original}{position}{to_aa} applied")
        return mutated

    @staticmethod
    def validate(sequence: str, position: int, from_aa: str) -> bool:
        idx = position - 1
        actual = sequence[idx]
        if actual != from_aa:
            log.warning(
                f"Position {position}: expected '{from_aa}' but sequence has '{actual}' "
                f"— check numbering offset"
            )
            return False
        return True

    @staticmethod
    def is_protected(position: int, protected_residues: set) -> bool:
        return position in protected_residues


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_mutation_string(mut_str: str) -> tuple:
    """Parse 'A316P' → ('A', 316, 'P')."""
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', mut_str.strip())
    if not m:
        raise ValueError(f"Cannot parse mutation '{mut_str}' — expected format like 'A316P'")
    return m.group(1), int(m.group(2)), m.group(3)


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

class RationalPipeline:

    def __init__(self, config: PipelineConfig, rounds: int, workdir: Path, template: Path = None):
        self.config = config
        self.rounds = rounds
        self.workdir = workdir
        workdir.mkdir(parents=True, exist_ok=True)

        self.boltz = BoltzRunner(config, workdir, template)
        self.claude = ClaudeReasoner(config)
        self.mutator = SequenceMutator()

        self.history: list = []

    def run(self):
        """
        Iterative discovery mode: Boltz-2 → Claude proposes mutation → validate → repeat.
        Requires Claude CLI and 'context' in config.
        """
        if not self.claude.is_available():
            log.error(
                "Iterative discovery mode requires the Claude CLI. "
                "Install it from https://claude.ai/code or use --scan_mutations "
                "to run a predefined mutation list without Claude."
            )
            sys.exit(1)

        if not self.config.context:
            log.error(
                "'context' is not set in the config YAML. "
                "Add a 'context' field describing the reaction mechanism and key residues."
            )
            sys.exit(1)

        log.info("=" * 60)
        log.info(f"{self.config.enzyme_name} × {self.config.substrate_name}  —  AI-Rational Pipeline")
        log.info(f"Rounds: {self.rounds}  |  Workdir: {self.workdir}")
        log.info("=" * 60)

        current_sequence = self.config.sequence
        mutations_applied: list = []
        carried_result: Optional[BoltzResult] = None

        for rnd in range(1, self.rounds + 1):
            log.info(f"\n{'─'*60}")
            log.info(f"ROUND {rnd}  |  Mutations so far: {mutations_applied or ['none']}")
            log.info(f"{'─'*60}")

            if carried_result is not None:
                boltz_result = carried_result
                carried_result = None
                log.info("  Reusing Boltz result carried from previous round")
            else:
                yaml_path = self.boltz.write_yaml(current_sequence, rnd)
                boltz_result = self.boltz.run(yaml_path, rnd)

            rationale = self.claude.reason(
                round_num=rnd,
                sequence=current_sequence,
                mutations_so_far=mutations_applied,
                boltz_result=boltz_result,
                history=self.history,
            )

            if rationale:
                log.info(f"\n  Claude reasoning:\n  {rationale.reasoning}")
                log.info(f"  Proposed: {rationale.proposed_mutation.get('mutation_string')}")
                log.info(f"  Expected: {rationale.expected_effect}")
                log.info(f"  Confidence: {rationale.confidence}")
                if rationale.warning:
                    log.warning(f"  Warning: {rationale.warning}")

            mutation_applied = None
            if rationale:
                prop = rationale.proposed_mutation
                position = prop.get("position")
                from_aa = prop.get("from_aa")
                to_aa = prop.get("to_aa")

                if position and to_aa:
                    if self.mutator.is_protected(position, self.config.protected_residues):
                        log.warning(f"  Proposed mutation at protected position {position} — skipping")
                    elif self.mutator.validate(current_sequence, position, from_aa):
                        new_sequence = self.mutator.apply(current_sequence, position, to_aa)
                        mut_str = prop.get("mutation_string", f"{from_aa}{position}{to_aa}")

                        log.info(f"  Validating {mut_str} with Boltz-2...")
                        yaml_mut = self.boltz.write_yaml(new_sequence, rnd, label=mut_str)
                        boltz_mutant = self.boltz.run(yaml_mut, rnd, label=mut_str)

                        delta = boltz_mutant.affinity_pred_value - boltz_result.affinity_pred_value
                        log.info(f"  Δaffinity = {delta:+.4f} ({'improved ✓' if delta < 0 else 'worse ✗'})")

                        if delta < 0 or rationale.confidence == "high":
                            current_sequence = new_sequence
                            mutations_applied.append(mut_str)
                            mutation_applied = mut_str
                            carried_result = boltz_mutant
                            log.info(f"  ✓ {mut_str} accepted")
                        else:
                            log.info(f"  {mut_str} not accepted (Δ={delta:+.4f}, confidence={rationale.confidence})")

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

    def run_mutation_scan(self, mutations: list, rollback_threshold: float = 0.3) -> list:
        """
        Scan a predefined ordered mutation list with Boltz-2.
        No Claude required. Claude commentary is added to step reports if available.

        For each mutation:
          - Apply to current best sequence
          - Run Boltz-2
          - Accept if Δaffinity ≤ rollback_threshold, else revert
        """
        log.info("=" * 60)
        log.info(f"{self.config.enzyme_name} × {self.config.substrate_name}  —  Mutation Scan")
        log.info(f"Mutations (in order): {mutations}")
        log.info(f"Rollback threshold  : Δ > +{rollback_threshold} log10(IC50) → skip")
        if self.claude.is_available() and self.config.context:
            log.info("  Claude CLI detected — will add structural commentary to step reports")
        log.info("=" * 60)

        current_sequence = self.config.sequence
        mutations_applied: list = []
        rows: list = []

        yaml_wt = self.boltz.write_yaml(current_sequence, 0, label="WT")
        current_result = self.boltz.run(yaml_wt, 0, label="WT")
        rows.append(self._scan_row(
            step=0, label="WT", mutations=[],
            result=current_result, delta=0.0, accepted=True, reason="baseline",
        ))

        for step, mut_str in enumerate(mutations, start=1):
            from_aa, position, to_aa = parse_mutation_string(mut_str)

            if SequenceMutator.is_protected(position, self.config.protected_residues):
                log.warning(f"  Step {step}: {mut_str} targets protected position {position} — skipped")
                rows.append(self._scan_row(
                    step=step, label=mut_str, mutations=list(mutations_applied),
                    result=current_result, delta=0.0, accepted=False, reason="protected residue",
                ))
                continue

            if not SequenceMutator.validate(current_sequence, position, from_aa):
                log.warning(f"  Step {step}: {mut_str} — residue mismatch, skipped")
                rows.append(self._scan_row(
                    step=step, label=mut_str, mutations=list(mutations_applied),
                    result=current_result, delta=0.0, accepted=False, reason="residue mismatch",
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
            log.info(f"  Δaffinity_pred_value = {delta:+.4f} (threshold = +{rollback_threshold})")

            if delta <= rollback_threshold:
                accepted = True
                reason = f"Δ={delta:+.4f} ≤ threshold"
                current_sequence = candidate_sequence
                current_result = mut_result
                mutations_applied.append(mut_str)
                log.info(f"  ✓ {mut_str} accepted")
            else:
                accepted = False
                reason = f"Δ={delta:+.4f} > threshold — rolling back"
                log.info(f"  ✗ {mut_str} skipped — {reason}")

            row = self._scan_row(
                step=step, label=cumulative_label, mutations=list(mutations_applied),
                result=mut_result, delta=delta, accepted=accepted, reason=reason,
            )
            rows.append(row)

            # Optional Claude commentary on this step
            rationale = self.claude.reason(
                round_num=step,
                sequence=current_sequence,
                mutations_so_far=list(mutations_applied),
                boltz_result=mut_result,
                history=[],
            )

            step_dir = self.workdir / f"scan_step{step}_{mut_str}"
            step_dir.mkdir(exist_ok=True)
            step_report = dict(row)
            if rationale:
                step_report["claude_reasoning"] = rationale.reasoning
                step_report["claude_expected_effect"] = rationale.expected_effect
                step_report["claude_confidence"] = rationale.confidence
                step_report["claude_alternatives"] = rationale.alternatives
                step_report["claude_warning"] = rationale.warning
            (step_dir / "step_report.json").write_text(json.dumps(step_report, indent=2))

        df = pd.DataFrame(rows)
        csv_path = self.workdir / "scan_summary.csv"
        df.to_csv(csv_path, index=False)

        seq_path = self.workdir / "scan_final_sequence.txt"
        label = "+".join(mutations_applied) or "WT"
        seq_path.write_text(f">{self.config.enzyme_name}_scan_{label}\n{current_sequence}\n")

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
    def _scan_row(step, label, mutations, result, delta, accepted, reason) -> dict:
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

        report: dict = {
            "round": s.round_num,
            "mutations_so_far": s.mutations_so_far,
            "mutation_applied_this_round": s.mutation_applied,
            "boltz": {
                "affinity_pred_value": s.boltz_result.affinity_pred_value,
                "affinity_probability": s.boltz_result.affinity_probability,
                "confidence": s.boltz_result.confidence_score,
                "contacts": s.boltz_result.contacts,
            },
        }
        if s.claude_rationale:
            report["claude_reasoning"] = s.claude_rationale.reasoning
            report["claude_proposal"] = s.claude_rationale.proposed_mutation
            report["expected_effect"] = s.claude_rationale.expected_effect
            report["claude_confidence"] = s.claude_rationale.confidence
            report["alternatives"] = s.claude_rationale.alternatives
            report["warning"] = s.claude_rationale.warning

        (d / "round_report.json").write_text(json.dumps(report, indent=2))
        log.info(f"  Saved round report → {d / 'round_report.json'}")

    def _save_final_report(self):
        rows = []
        for s in self.history:
            row = {
                "round": s.round_num,
                "mutation_applied": s.mutation_applied or "none",
                "cumulative_mutations": "+".join(s.mutations_so_far) or "WT",
                "affinity_pred_value": round(s.boltz_result.affinity_pred_value, 4),
                "binder_probability": round(s.boltz_result.affinity_probability, 4),
                "confidence": round(s.boltz_result.confidence_score, 4),
                "n_contacts": len(s.boltz_result.contacts),
                "contacts": ", ".join(s.boltz_result.contacts),
            }
            if s.claude_rationale:
                row["claude_confidence"] = s.claude_rationale.confidence
            rows.append(row)

        df = pd.DataFrame(rows)
        csv_path = self.workdir / "pipeline_summary.csv"
        df.to_csv(csv_path, index=False)

        log.info("\n" + "=" * 60)
        log.info("EVOLUTION TRAJECTORY")
        log.info("=" * 60)
        log.info(df[["round","mutation_applied","cumulative_mutations",
                      "affinity_pred_value","binder_probability"]].to_string(index=False))

        if self.history:
            final_seq = self.history[-1].sequence
            final_muts = self.history[-1].mutations_so_far
            seq_path = self.workdir / "final_sequence.txt"
            seq_path.write_text(
                f">{self.config.enzyme_name}_optimised_{'_'.join(final_muts) or 'WT'}\n{final_seq}\n"
            )
            log.info(f"\nFinal sequence → {seq_path}")
            log.info(f"Mutations: {final_muts or 'none'}")

        log.info(f"Summary CSV → {csv_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Enzyme × substrate AI-rational iterative mutation pipeline"
    )
    p.add_argument(
        "--config", type=Path, required=True,
        help="Path to experiment config YAML (see configs/pada1_nnbt.yaml for example)"
    )
    p.add_argument("--rounds", type=int, default=4,
        help="Number of rounds (iterative discovery mode only)")
    p.add_argument("--workdir", type=Path, default=Path("results"))
    p.add_argument(
        "--template", type=Path, default=None,
        help="CIF/PDB file for backbone template constraint"
    )
    p.add_argument(
        "--scan_mutations", type=str, default=None,
        help=(
            "Comma-separated ordered mutations to scan (e.g. 'A316P,F191L,A73T'). "
            "Activates scan mode — no Claude required. "
            "Each mutation is tested in order; if it worsens affinity beyond "
            "--rollback_threshold it is skipped."
        ),
    )
    p.add_argument(
        "--rollback_threshold", type=float, default=0.3,
        help="Max allowed Δaffinity for acceptance in scan mode (default: 0.3)"
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    config = load_config(args.config)

    pipeline = RationalPipeline(
        config=config,
        rounds=args.rounds,
        workdir=args.workdir,
        template=args.template,
    )

    if args.scan_mutations:
        mutations = [m.strip() for m in args.scan_mutations.split(",") if m.strip()]
        pipeline.run_mutation_scan(mutations, rollback_threshold=args.rollback_threshold)
    else:
        pipeline.run()
