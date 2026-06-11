
# MutateX - FoldX

## Paper overview & explanation 

A pipeline that automates in silico saturation mutagenesis by systematically substituting every residue in a protein to all other amino acids and computing the free energy change (ΔΔG) for each. MutateX computes this using FoldX's empirical energy function and FoldX models contributions from main-chain and side-chain interactions.

ΔΔG = free energy of the mutant − free energy of the wild-type
Positive value = the mutation is destabilizing (the mutant is less stable) 
Negative value = it's stabilizing (the mutant is more stable)

**FoldX** works in three steps for each mutation: 
(1) Repair — energy-minimizes the input structure to remove clashes and optimize side chains
(2) BuildModel — generates the mutant and a matching wild-type model by swapping the side chain 
(3) Energy calculation — calculates ΔΔG as the difference. For binding affinity, it also runs AnalyseComplex

Key limitation: FoldX does not move the protein backbone, but only rearranges side chains locally. This means it can overestimate the cost of inserting a large residue where a small one was, because there's no room to adapt + it also means long-range allosteric effects aren't captured.

Saturation = every single residue in the protein is mutated to every other amino acid, resulting in a complete "mutational landscape" -> a map of how sensitive each position is to change. This can reveal which residues are structurally critical (high average ΔΔG) which are tolerant to substitution (ΔΔG ≈ 0), and which rare positions where a mutation actually improves stability.

!! empirical !! — the energy function was trained by fitting coefficients to a database of experimentally measured ΔΔG values from single-point mutants. This means the values it produces are not physical free energies in an absolute sense, but they are well-calibrated to match experimental trends.

The three-step mutation protocol:

1) Repair (RepairPDB)

Before doing anything, FoldX energy-minimizes the input structure. This removes steric clashes that exist in the raw PDB file and optimizes rotamers (side-chain conformations). Without this, small imperfections in the crystal structure would produce noisy ΔΔG values. MutateX runs Repair once per input model before any mutations are attempted.

2) Build model (BuildModel)

FoldX swaps the specified residue's side chain for the new amino acid and selects the best rotamer from a library. Critically, the protein backbone is kept completely fixed — only side chains move. A matching wild-type model is also generated (same backbone relaxation, original residue) so the comparison is fair. For homo-multimers, MutateX makes the same mutation on all identical chains simultaneously.

3) Energy calculation (and optionally AnalyseComplex)

FoldX evaluates its empirical energy function on both the mutant and matched wild-type models and subtracts: ΔΔG = ΔG(mutant) − ΔG(WT). For protein complexes, AnalyseComplex additionally computes the free energy of interaction between chains, giving a binding ΔΔG on top of the folding ΔΔG.

**The energy function**

The FoldX energy function sums several physical and statistical terms. Each is weighted by a coefficient fitted to the experimental database:

ΔG = w₁·ΔGvdW + w₂·ΔGsolvH + w₃·ΔGsolvP + w₄·ΔGwb + w₅·ΔGhbond + w₆·ΔGel + w₇·ΔGkon + w₈·ΔGmc + w₉·ΔGsc

ΔGvdW = Van der Waals interactions — steric packing between atoms

ΔGsolvH / ΔGsolvP = Solvation energies for hydrophobic and polar groups — burial of hydrophobic surface drives folding

ΔGhbond = Hydrogen bond energy — counts and quality of H-bonds made

ΔGel = Electrostatic interactions — charges, salt bridges

ΔGwb = Water bridges — structural water molecules mediating interactions

ΔGkon = Entropic cost — conformational entropy change on mutation

ΔGmc / ΔGsc = Main-chain and side-chain entropy terms



LIMITS: What FoldX cannot predict

*Allosteric effects* — FoldX only samples local side-chain space. A mutation 20 Å away from the active site that shifts domain dynamics will not be detected.
*Epistasis — each mutation is treated independently. Cooperative effects where two distant mutations interact are invisible (e.g. Q498R + N501Y in the COVID spike work together but FoldX can't model that).
*Pi-stacking* — the FoldX force field does not properly represent aromatic ring stacking interactions. This is why N501Y in Omicron (stabilizing via pi-stacking with Y41 on ACE2) was incorrectly predicted as destabilizing.
*Backbone conformational changes* — if a mutation causes a loop to remodel or a helix to shift, FoldX will miss this entirely. The backbone is frozen.
*Functional effects beyond stability* — a mutation can disrupt catalysis, ligand binding, or a protein-protein interaction without destabilizing the fold. FoldX only sees stability and direct binding energetics.





'''
this is code
'''