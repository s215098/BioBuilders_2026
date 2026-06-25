# How to Use Claude + Claude Code to Build OT-2 Protocols

*iGEM 2026 — internal working guide. Commit this to the repo; it's for the whole team.*

This explains how we turn an assay (e.g. the NNBT colorimetric screen) into a **validated, LP-compliant OT-2 protocol** using two different Claude surfaces, and how we test it from dry simulation through to a real run without wasting reagents or robot time.

---

## 1. The mental model: two Claudes, two jobs

The single most important idea: **the chat does the thinking, Claude Code does the building.** Keep them in their lanes.

| | **Claude chat** (claude.ai project) | **Claude Code** (in the repo) |
|---|---|---|
| Lives in | The project, with the brief + paper summaries + this guide as context | Your GitHub repo, with the LP convention + existing protocols + `CLAUDE.md` |
| Job | *Design & review*: assay logic, automate-vs-manual, deck layout, parameters, acceptance criteria, scientific sanity | *Implement & validate*: write the `.ipynb`/`.py`, run the validation gates, fix errors, commit, open PRs |
| Decides | The science (volumes, conditions, controls — sourced from papers) | Nothing about the science; it executes the spec |
| Output | A **protocol design spec** (a structured handoff) | A simulated, syntax-clean, committed protocol |
| Touches the repo? | No | Yes |

The reason for the split: the chat has the biology context and won't hallucinate Opentrons API calls into your repo; Claude Code has the repo context and the validation tooling but should **never invent assay parameters**. If Code needs a number that isn't in the spec, that's a signal to come back to the chat, not to guess.

---

## 2. One-time setup (do this once per repo)

1. **Add a `CLAUDE.md` at the repo root.** Claude Code reads it automatically on every session — this is the highest-leverage thing you can do. A starter is in §7 below. It pins the ID scheme, the validation gates, and the hardware so Code is compliant by default.
2. **Keep one "golden" protocol** in the repo as a template. Once `W-O-SC-01` passes review, it becomes the reference every future protocol is cloned from. Point `CLAUDE.md` at it.
3. **Install the simulator** in the repo's environment: `pip install opentrons`. This gives you `opentrons_simulate`, which is a required gate in our LP convention.
4. **Branch-per-workflow.** One branch per workflow ID (`feat/W-O-SC-01`), one PR, reviewed by a second person → that reviewer is the "Approved by" in the SOP header.

---

## 3. The per-protocol loop

This is the cadence for every protocol. Steps marked **[chat]** happen here in the project; **[code]** happen in Claude Code; **[lab]** happen at the bench.

1. **[chat] Design session.** We work out the assay → automate/manual split → deck layout → parameters → acceptance criteria. Output: a **protocol design spec** (see §4).
2. **[code] Scaffold.** In the repo: paste the spec, point Code at the LP convention + the golden example, ask it to scaffold the notebook for workflow ID `W-O-…`.
3. **[code] Validation gates.** Code runs, in order: `py_compile` (syntax) → `opentrons_simulate` (deck conflicts, tip math, volumes, collisions). It iterates until both pass. **No green gates, no lab.**
4. **[chat] Scientific review.** Bring the simulation summary (and any "Code asked me for a number I didn't have" moments) back here. We check the science, not the syntax.
5. **[code] Commit + PR.** Branch, commit, PR. Second teammate reviews and approves.
6. **[lab] Dry-to-wet testing ladder** (see §5).
7. **[chat/code] Capture deviations** into the Run Record + change log, and improve the protocol.

The loop is deliberately short and repeatable. Most protocols will go around it 2–3 times before a wet run.

---

## 4. The handoff artifact: the "protocol design spec"

This is what the chat produces and Claude Code consumes. Keeping it explicit is what stops Code from guessing. A good spec contains:

- **Workflow ID + name** (e.g. `W-O-SC-01`, NNBT colorimetric screen)
- **Unit operations**, in order, each as a one-line goal (these become the table-of-contents)
- **Parameters**, fully named with units (e.g. `SUPERNATANT_VOL_UL = 100`) — the single source of truth, no magic numbers
- **Deck layout**: which labware in which slot, which pipette does what, reservoir/reagent map
- **Liquid map**: source wells/reservoir → destination wells, per step
- **Acceptance criteria** per unit op (what must be true to proceed)
- **Controls + plate layout** (where the parental standard / reference enzyme / blanks go)
- **Off-deck handoffs** (e.g. "read at 530 nm on the Synergy HT" → a `U-C-` plate-reader unit op, not robot code)

If a number isn't in the spec, Code shouldn't be inventing it — that's the guardrail.

---

## 5. The testing ladder (dry → wet)

Never jump from code to a full reagent run. Climb the ladder; each rung is cheaper than the one above it and catches a different class of error.

1. **`py_compile`** — syntax only. Seconds, no robot. *Catches:* typos, bad indentation.
2. **`opentrons_simulate`** — full protocol simulation, no robot. *Catches:* deck collisions, running out of tips, volume/labware mismatches, impossible transfers. **Mandatory gate.**
3. **Deck-layout review** — eyeball the slot map and reagent volumes against the spec. *Catches:* reservoir fill too low (dead volume), wrong plate in a slot.
4. **Labware Position Check (LPC)** — on the real OT-2, in the Opentrons App, *no liquid*. Calibrate and save labware offsets. *Catches:* tips crashing into plate edges, mis-seated labware. Export the offsets and keep them with the protocol.
5. **Water / dye dry run** — run the *whole* protocol with water (add food dye for visibility). Optionally weigh plates before/after for a gravimetric volume check. *Catches:* wrong wells, spills, aspiration-height problems, mixing issues, timing. This is the rung most people skip and most regret skipping.
6. **Reagent test run — one plate, with controls.** Real reagents, full controls (parental standard, reference enzyme, blanks). Compare against the manual benchmark. For the NNBT screen, this is also where you check your CV against Dolz's 19% manual baseline — automation should *improve* it.
7. **Full production run.**

Rules of thumb: rungs 1–3 are free and fast — there is no excuse to skip them. Rungs 4–5 cost only water and an hour of robot time and save you reagents and failed screens. Only spend real reagents (rung 6+) on a protocol that passed water.

> **NNBT-specific cautions for the wet rungs:** the Purpald reagent is 50 mM in **2N NaOH** — caustic; dedicate tips and rinse the deck after. Watch the NNBT-in-acetonitrile addition order (precipitation risk per Dolz). Use the temperature module to keep the H₂O₂/NNBT mix cold during the run, and the heater-shaker for both the 30 °C/5 min incubation and the 10 min Purpald shake.

---

## 6. Guardrails (read before every wet run)

- **The robot faithfully executes a wrong protocol.** A human reviews and approves every protocol before the first wet run. No exceptions.
- **Code builds, chat decides the science.** If Claude Code proposes an assay parameter, stop — that number belongs in the spec, sourced from a paper or a design decision, not from Code.
- **Pin the Opentrons API level** in the protocol metadata so the simulator and the robot agree.
- **Single source of truth for numbers.** Everything lives in the Parameters section; if you change a volume, change it there, re-simulate, re-commit.
- **Every run writes a Run Record** (operator, date, pre/during/deviation/post notes) — it's in the LP template and it's what feeds the change log.

---

## 7. Starter `CLAUDE.md` (drop in repo root)

```markdown
# OT-2 Protocol Repo — Claude Code instructions

## What this repo is
LP-compliant Opentrons OT-2 protocols for the iGEM 2026 epoxy-degradation project.

## Hard rules
- Follow the literate programming convention in `docs/01_How_to_Write_a_Literate_Programming_Protocol.docx`.
- One notebook = one workflow. ID format: W-O-TT-XX (TT: EN=enzymatic, TR=transformation,
  PU=purification, SC=screening). Unit ops: U-O-XX (Opentrons), U-C-XX (Clariostar/reader).
- All numbers live in ONE Parameters section, named with units (e.g. PURPALD_VOL_UL). No magic numbers.
- Each unit op section: Goal / Inputs / Outputs / Acceptance criteria / Record / Deviation handling / Notes.
- Separate markdown narrative from code cells.
- Do NOT invent assay parameters. If a value is missing, stop and ask — it comes from the design spec.

## Validation gates (must all pass before PR)
- `python -m py_compile <generated_protocol>.py`
- `opentrons_simulate <generated_protocol>.py`
- Pin the API level in protocol metadata.

## Hardware available
- Pipettes: P300 8-channel, P300 single, P20 8-channel. (No P20 single-channel — keep sub-20 µL columnar.)
- Modules: heater-shaker, temperature, thermocycler (PCR), HEPA.
- No on-deck centrifuge, colony picker, or electroporator.

## Golden example
Clone the structure of `protocols/W-O-SC-01_nnbt_screen/` for any new protocol.
```

---

## 8. Worked example: building `W-O-SC-01` (NNBT screen)

1. **[chat]** We finalise the spec: unit ops = supernatant transfer → add NNBT mix (30 °C, 5 min, shake) → add Purpald (shake 10 min) → handoff to reader. Parameters, deck map, plate layout with the parental standard column and reference-enzyme well, acceptance criteria.
2. **[code]** In the repo, on branch `feat/W-O-SC-01`: paste spec, point Code at `CLAUDE.md` + the LP doc, ask it to scaffold the notebook that emits `W-O-SC-01_nnbt_screen.py`.
3. **[code]** Code runs `py_compile` then `opentrons_simulate`, fixes tip/volume errors, gets both green.
4. **[chat]** We review the simulation summary — tip counts, reservoir fills (dead volume!), well mapping — and the plate layout for the science.
5. **[code]** Commit, PR, teammate approves (→ "Approved by").
6. **[lab]** LPC (no liquid) → water+dye dry run (check wells, no spills, timing) → one-plate reagent test with controls, read at 530 nm, compare CV to Dolz baseline → full run.
7. **[chat/code]** Log deviations, update change log.

---

## Quick-reference checklist

- [ ] Design spec produced in chat (parameters, deck, layout, acceptance criteria, controls)
- [ ] `CLAUDE.md` and golden example in place
- [ ] `py_compile` passes
- [ ] `opentrons_simulate` passes
- [ ] Deck layout + reagent volumes reviewed
- [ ] PR reviewed and approved by a second person
- [ ] LPC done on robot (offsets saved)
- [ ] Water/dye dry run clean
- [ ] One-plate reagent test with controls passes
- [ ] Run Record + change log updated
