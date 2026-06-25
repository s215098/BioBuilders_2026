"""
=======================================================================================
OT2 Protocol: UPO_GGA_Transformation (W-O-CL-01)
=======================================================================================
 Doc No:        W-O-CL-01
 Version:       v4.1
 Date:          2026-06-22
 Owner:         Lukas 
 Workflow ID:   W-O-CL-01
 Unit ops:      U-O-01, U-O-02, U-O-03, U-O-04, U-T-01, U-M-01, U-O-05, U-T-02, U-M-02
 Output:        N chemically-transformed E. coli cultures + 1 no-insert negative
                 control (post heat-shock, pre-recovery-plating), ready for manual
                 plating onto individual agar plates.

   T = Thermocycler unit operation
   M = Manual (human) unit operation - protocol pauses for an off-robot step


 SCOPE / DESIGN RATIONALE (carried over):
 - Backbone is a single standard vector, identical for all reactions. Each
   UPO is supplied as ONE gBlock (promoter + terminator + restriction sites +
   gene, all on one fragment) -> a 2-fragment GG reaction per construct
   (backbone + 1 gBlock).
 - Plating is fully manual (U-M-02): no 6-well (or other) agar plate is loaded.
 - Tube rack capacity note: a single 24-tube rack holds backbone + master mix
   (1 tube each) + up to 22 gBlocks. If NUM_UPOS grows beyond 22, a second
   tube rack (and source-lookup logic) would be needed -- not handled here.
=======================================================================================
 STEPS TABLE
=======================================================================================
 U-O-01   Distribute backbone from its tube (P20, single-nozzle)
 U-O-02   Distribute master mix from its tube (P20, single-nozzle)
 U-O-03   Distribute water (P20, full 8-channel, reservoir)
 U-O-04   Add the 20 gBlocks + negative-control water, one tube/well at a time
 U-T-01   Golden Gate cycling program (thermocycler)
 U-M-01   PAUSE - operator loads cold competent cells into the reservoir
 U-O-05   Distribute + gently mix competent cells (P300, full 8-channel)
 U-T-02   Heat-shock + recovery program (thermocycler)
 U-M-02   PAUSE - protocol ends; operator manually plates all wells
=======================================================================================
"""

import math
from opentrons import protocol_api
from opentrons.protocol_api import SINGLE

metadata = {
    'apiLevel': '2.20',
    'protocolName': 'W-O-CL-01 - UPO Golden Gate Assembly + Heat-Shock Transformation',
    'description': (
        'GG assembly of a standard pPIC backbone with N gBlocks plus a no-insert '
        'negative control -- all pipetted automatically via single-nozzle partial '
        'tip pickup from individual tubes -- followed by competent-cell mixing and '
        'heat-shock transformation. Plating is manual.'
    ),
}

# =======================================================================================
# PARAMETERS  (single source of truth - no magic numbers anywhere else in this file)
# =======================================================================================

NUM_UPOS = 20                       # distinct UPO gBlock constructs -- change this and
                                     # REACTION_COLUMNS / GBLOCK_TUBE_MAP / NEGATIVE_
                                     # CONTROL_WELL all recompute automatically below
                                     # (tube rack holds up to 22 gBlocks -- see note above)
WELLS_PER_COLUMN = 8
NUM_NEEDED_WELLS = NUM_UPOS + 1      # +1 for the no-insert negative control
NUM_REACTION_COLUMNS = math.ceil(NUM_NEEDED_WELLS / WELLS_PER_COLUMN) #math.ceil rounds up to the next whole column, since a partial column is still a full column on the plate
REACTION_COLUMNS = [f'A{i + 1}' for i in range(NUM_REACTION_COLUMNS)] #for example, if NUM_REACTION_COLUMNS = 3, this will generate ['A1', 'A2', 'A3']
ROW_LETTERS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Per-well GG reaction recipe (uL) -- robot adds all 20 uL now, no manual addition step -> can be changed by user
VOL_MASTER_MIX_UL = 10.0            # 2X Golden Gate enzyme + buffer master mix
VOL_WATER_UL = 8.0
VOL_BACKBONE_UL = 1.0               # standard pPIC backbone, normalized stock
VOL_GBLOCK_UL = 1.0
GG_REACTION_VOLUME_UL = VOL_MASTER_MIX_UL + VOL_WATER_UL + VOL_BACKBONE_UL + VOL_GBLOCK_UL  # 20 uL

# Heat-shock transformation recipe
VOL_COMPETENT_CELLS_UL = 50.0
HS_REACTION_VOLUME_UL = GG_REACTION_VOLUME_UL + VOL_COMPETENT_CELLS_UL  # 70 uL
MIX_AFTER_REPS = 3
MIX_AFTER_VOL_UL = 30.0
MIX_FLOW_RATE_UL_PER_S = 40.0       # slowed from the P300's default 94 uL/s -- gentler
                                     # on shear-sensitive competent cells

# Golden Gate thermocycler profile (max sample temp in this profile: 60 degC)
GG_CYCLE_TEMP_HIGH_C = 37
GG_CYCLE_TIME_HIGH_S = 120
GG_CYCLE_TEMP_LOW_C = 16
GG_CYCLE_TIME_LOW_S = 300
GG_CYCLE_REPETITIONS = 25
GG_DEACTIVATION_TEMP_C = 60
GG_DEACTIVATION_TIME_S = 300
GG_HOLD_TEMP_C = 4
GG_LID_TEMP_C = 70                  # ~10 degC above this profile's 60 degC max, not the
                                     # PCR-style 105 degC the original script copied over

# Heat-shock thermocycler profile (max sample temp in this profile: 42 degC)
HS_COLD_SOAK_TEMP_C = 4
HS_COLD_SOAK_TIME_S = 600
HS_SHOCK_TEMP_C = 42
HS_SHOCK_TIME_S = 90
HS_COLD_RECOVER_TEMP_C = 4
HS_COLD_RECOVER_TIME_S = 120
HS_RECOVERY_TEMP_C = 37
HS_RECOVERY_TIME_S = 3600
HS_LID_TEMP_C = 47                  # ~5 degC above this profile's 42 degC max
HS_FINAL_HOLD_TEMP_C = 4            # block holds here through the U-M-02 plating pause,
                                     # keeping cells cold until plating instead of
                                     # drifting to ambient after the run ends

# Single-nozzle dispensing chunk size: how much to aspirate per refill when
# distributing a uniform reagent (backbone, master mix) with one reused tip.
# Kept comfortably under the P20's 20 uL max to leave headroom for the disposal/
# overage the API adds on each aspirate.
SINGLE_NOZZLE_CHUNK_UL = 18.0

# Tube rack (opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap, slot 5) has
# 4 ROWS (A-D), unlike the reaction plate's 8 -- a separate row count is needed here
# to avoid generating nonexistent well names like 'E1'.
TUBERACK_ROWS_PER_COLUMN = 4
TUBERACK_BACKBONE_WELL = 'A6'
TUBERACK_MASTERMIX_WELL = 'B6'

GBLOCK_TUBE_MAP = {}                 # tube-rack well -> gBlock name
for _i in range(NUM_UPOS):
    _col = _i // TUBERACK_ROWS_PER_COLUMN
    _row = ROW_LETTERS[_i % TUBERACK_ROWS_PER_COLUMN]
    GBLOCK_TUBE_MAP[f'{_row}{_col + 1}'] = f'UPO_{_i + 1:02d}' 

# gBlock -> reaction-plate well map (reaction-side; independent of tube layout above)
GBLOCK_WELL_MAP = {}
for _i in range(NUM_UPOS):
    _col = _i // WELLS_PER_COLUMN
    _row = ROW_LETTERS[_i % WELLS_PER_COLUMN]
    GBLOCK_WELL_MAP[f'{_row}{_col + 1}'] = f'UPO_{_i + 1:02d}'

_neg_col = NUM_UPOS // WELLS_PER_COLUMN
_neg_row = ROW_LETTERS[NUM_UPOS % WELLS_PER_COLUMN]
NEGATIVE_CONTROL_WELL = f'{_neg_row}{_neg_col + 1}'

# Reservoir well assignments (nest_12_reservoir_15ml) -- water (bulk, also used for
# the single negative-control well) and competent cells (bulk, loaded mid-run)
RESERVOIR_WATER_WELL = 'A1'
RESERVOIR_CELLS_WELL = 'A2'

# Deck slots. Slots 4 and 5 are used for the
# tip rack and tube rack specifically because Opentrons' partial-tip-pickup deck
# extents documentation confirms these two slots have no edge-of-deck restrictions
# for single-nozzle access (front-row slots 1-3 do; using slot 2 caused a real
# "outside of robot bounds" error during validation).
SLOT_TUBERACK = '5'
SLOT_RESERVOIR = '3'
SLOT_TIPRACK_20 = '4'
SLOT_TIPRACK_300 = '6'


def _ordered_reaction_well_names(num_wells):
    """First `num_wells` well names of the reaction plate, in column-then-row order
    (A1, B1, ... H1, A2, B2, ...) -- matches GBLOCK_WELL_MAP / NEGATIVE_CONTROL_WELL."""
    names = []
    for col in range(NUM_REACTION_COLUMNS):
        for row in ROW_LETTERS:
            names.append(f'{row}{col + 1}')
            if len(names) == num_wells:
                return names
    return names


# =======================================================================================
def run(protocol: protocol_api.ProtocolContext):

    # -----------------------------------------------------------------------------
    # DECK + INSTRUMENT SETUP
    # -----------------------------------------------------------------------------
    tr_20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    tr_300 = protocol.load_labware('opentrons_96_tiprack_300ul', SLOT_TIPRACK_300)

    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr_20])
    p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tr_300])

    # P20 stays in single-nozzle mode for its entire role in this protocol (backbone,
    # master mix, gBlocks, negative control all use it) -- configured once here.
    p20.configure_nozzle_layout(style=SINGLE, start='A1', tip_racks=[tr_20])

    tc_mod = protocol.load_module('thermocycler module')
    reaction_plate = tc_mod.load_labware('biorad_96_wellplate_200ul_pcr') #tbd
    tc_mod.open_lid()
    tc_mod.set_block_temperature(GG_HOLD_TEMP_C)

    tuberack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', #TBD
                                      SLOT_TUBERACK)
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', SLOT_RESERVOIR)

    reaction_well_names = _ordered_reaction_well_names(NUM_NEEDED_WELLS)
    reaction_wells = [reaction_plate[w] for w in reaction_well_names]

    # Column-anchor wells (top 'A'-row well of each used column) for full-8-channel
    # steps. A multichannel pipette addresses whole columns via these anchors, NOT the
    # individual wells in reaction_wells. REACTION_COLUMNS == ['A1', 'A2', 'A3'].
    # Note: the last column is addressed in full, so its unused wells (F3/G3/H3) also
    # receive water and cells -- harmless, but see competent-cell volume note in header.
    reaction_columns = [reaction_plate[w] for w in REACTION_COLUMNS]

    def distribute_single_nozzle(volume_per_well, source_well, dest_wells):
        """Distribute a uniform reagent to many wells with ONE reused tip, re-aspirating
        from source_well whenever the tip runs low. Safe to reuse the tip throughout:
        every destination receives the same reagent, so there's no cross-contamination
        concern the way there would be with different reagents."""
        p20.pick_up_tip()
        remaining = 0.0
        for well in dest_wells:
            if remaining < volume_per_well:
                p20.aspirate(SINGLE_NOZZLE_CHUNK_UL - remaining, source_well)
                remaining = SINGLE_NOZZLE_CHUNK_UL
            p20.dispense(volume_per_well, well)
            remaining -= volume_per_well
        p20.drop_tip()

    # ===============================================================================
    # U-O-01  Distribute backbone (P20, single-nozzle, one tube)
    # ===============================================================================
    distribute_single_nozzle(VOL_BACKBONE_UL, tuberack[TUBERACK_BACKBONE_WELL], reaction_wells)

    # ===============================================================================
    # U-O-02  Distribute master mix (P20, single-nozzle, one tube)
    # ===============================================================================
    distribute_single_nozzle(VOL_MASTER_MIX_UL, tuberack[TUBERACK_MASTERMIX_WELL], reaction_wells)

    # ===============================================================================
    # U-O-03  Distribute water (P20, full 8-channel, reservoir)
    #   This one stays full-8-channel: it's the same reagent across whole reaction
    #   columns sourced from a true reservoir channel, which is exactly what full
    #   8-channel access is efficient and safe for.
    # ===============================================================================
    # NOTE: the same physical P20 is reconfigured between single-nozzle and
    # full-8-channel layouts at the points where each is needed -- there's only one
    # P20 pipette attached, so this is a nozzle-layout switch, not a second pipette.
    p20.configure_nozzle_layout(style=protocol_api.ALL, tip_racks=[tr_20])
    p20.distribute(VOL_WATER_UL, reservoir[RESERVOIR_WATER_WELL], reaction_columns)
    p20.configure_nozzle_layout(style=SINGLE, start='A1', tip_racks=[tr_20])

    # ===============================================================================
    # U-O-04  Add the 20 gBlocks + negative-control water (P20, single-nozzle,
    #          one tube/well + one fresh tip per item)
    #   Goal:        1 uL of each gBlock into its mapped reaction well; 1 uL of plain
    #                water (not gBlock) into NEGATIVE_CONTROL_WELL.
    #   Acceptance:  All NUM_UPOS gBlocks delivered to their correct wells; exactly
    #                one well dosed with water instead, as the negative control.
    #   Note:        a fresh tip per item, unlike backbone/master mix above -- these
    #                are different reagents, and cross-contamination between
    #                different UPO constructs would compromise the comparison this
    #                whole experiment exists to make.
    # ===============================================================================
    for tube_well, gblock_name in GBLOCK_TUBE_MAP.items():
        dest_well = [w for w, n in GBLOCK_WELL_MAP.items() if n == gblock_name][0]
        p20.pick_up_tip()
        p20.aspirate(VOL_GBLOCK_UL, tuberack[tube_well])
        p20.dispense(VOL_GBLOCK_UL, reaction_plate[dest_well])
        p20.drop_tip()

    p20.pick_up_tip()
    p20.aspirate(VOL_GBLOCK_UL, reservoir[RESERVOIR_WATER_WELL])
    p20.dispense(VOL_GBLOCK_UL, reaction_plate[NEGATIVE_CONTROL_WELL])
    p20.drop_tip()

    # ===============================================================================
    # U-T-01  Golden Gate cycling (thermocycler)
    #   Lid temperature is set just above this profile's own max sample temperature
    #   (60 degC deactivation step), not the 105 degC PCR-style default.
    # ===============================================================================
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(GG_LID_TEMP_C)
    gg_profile = [
        {'temperature': GG_CYCLE_TEMP_HIGH_C, 'hold_time_seconds': GG_CYCLE_TIME_HIGH_S},
        {'temperature': GG_CYCLE_TEMP_LOW_C, 'hold_time_seconds': GG_CYCLE_TIME_LOW_S},
    ]
    tc_mod.execute_profile(steps=gg_profile, repetitions=GG_CYCLE_REPETITIONS,
                            block_max_volume=GG_REACTION_VOLUME_UL)
    tc_mod.execute_profile(
        steps=[{'temperature': GG_DEACTIVATION_TEMP_C, 'hold_time_seconds': GG_DEACTIVATION_TIME_S}],
        repetitions=1, block_max_volume=GG_REACTION_VOLUME_UL)
    tc_mod.set_block_temperature(GG_HOLD_TEMP_C)
    tc_mod.open_lid()

    # ===============================================================================
    # U-M-01  MANUAL - operator loads cold competent cells
    #   Reason: loading cells just before use (rather than at protocol start) keeps
    #           them cold/fresh through the ~1 h GG cycling step above.
    # ===============================================================================
    protocol.pause(
        'U-M-01 MANUAL STEP: place freshly-thawed, on-ice competent cells into '
        f'reservoir well {RESERVOIR_CELLS_WELL}, then resume.'
    )

    # ===============================================================================
    # U-O-05  Distribute + gently mix competent cells (P300, full 8-channel)
    #   Note: flow rate is temporarily slowed for the mix_after step, since
    #         competent cells are shear-sensitive; restored afterward.
    #   Targets column anchors (reaction_columns), not individual wells. new_tip is
    #   'always' so each column gets fresh tips: mix_after stirs DNA-containing wells,
    #   and reusing tips between columns would carry one construct into the next.
    # ===============================================================================
    default_aspirate_rate = p300.flow_rate.aspirate
    default_dispense_rate = p300.flow_rate.dispense
    p300.flow_rate.aspirate = MIX_FLOW_RATE_UL_PER_S
    p300.flow_rate.dispense = MIX_FLOW_RATE_UL_PER_S
    p300.transfer(VOL_COMPETENT_CELLS_UL, reservoir[RESERVOIR_CELLS_WELL], reaction_columns,
                  new_tip='always', mix_after=(MIX_AFTER_REPS, MIX_AFTER_VOL_UL))
    p300.flow_rate.aspirate = default_aspirate_rate
    p300.flow_rate.dispense = default_dispense_rate

    # ===============================================================================
    # U-T-02  Heat-shock + recovery (thermocycler)
    #   Lid temperature is set just above this profile's own max sample temperature
    #   (42 degC shock step), not inherited from the GG profile's 70 degC setting.
    # ===============================================================================
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(HS_LID_TEMP_C)
    hs_profile = [
        {'temperature': HS_COLD_SOAK_TEMP_C, 'hold_time_seconds': HS_COLD_SOAK_TIME_S},
        {'temperature': HS_SHOCK_TEMP_C, 'hold_time_seconds': HS_SHOCK_TIME_S},
        {'temperature': HS_COLD_RECOVER_TEMP_C, 'hold_time_seconds': HS_COLD_RECOVER_TIME_S},
        {'temperature': HS_RECOVERY_TEMP_C, 'hold_time_seconds': HS_RECOVERY_TIME_S},
    ]
    tc_mod.execute_profile(steps=hs_profile, repetitions=1, block_max_volume=HS_REACTION_VOLUME_UL)
    tc_mod.open_lid()
    # Hold the block cold through the U-M-02 pause (lid already open so the operator
    # can lift the plate out). Cells stay at HS_FINAL_HOLD_TEMP_C until plating rather
    # than drifting to ambient as they would after deactivate(). Brief condensation on
    # the open cold block during the plating window is acceptable. deactivate() is
    # deferred until after the operator finishes plating (below).
    tc_mod.set_block_temperature(HS_FINAL_HOLD_TEMP_C)

    # ===============================================================================
    # U-M-02  MANUAL - operator plates all wells (20 UPOs + 1 negative control)
    # ===============================================================================
    protocol.pause(
        f'PROTOCOL COMPLETE (U-M-02): manually plate all {NUM_UPOS} transformations '
        f'plus the negative control (well {NEGATIVE_CONTROL_WELL}) onto individual '
        'agar plates. The negative control should show few/no colonies -- treat '
        'real-well results with caution if it doesn\'t.'
    )

    # Operator has plated; release the cold hold and end the run.
    tc_mod.deactivate()
