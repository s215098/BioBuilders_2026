# Investigating Laccases from Geobacillus as potential enzymes for breakdown of BADGE

From this paper: https://doi.org/10.1007/s11356-024-34095-y (available in Google Drive under Research Articles) we know that laccases from Geobacillus thermoparafinivorans and Geobacillus stearothermoparafinivorns can degrade BADGE (Bisphenol A glycidyl ether). See Table 1 in the paper or figure below (last 2 entries). However Geobacillus stearothermoparafinivorns is mentioned in the title of the paper as Geobacillus stearothermophilus.

![](figures/laccase-degradation-bpa-vs-badge.png)

Next step is to search for these enzymes. Authors mention that they used NCBI, we are gonna start with UniProt

1. Head to: https://www.uniprot.org/

2. Change from UniProt to UniPark

3. By filling these fields in the "Advanced" Search Bar we see that information is only available for Geobacillus Stearothermophilus:

![](figures/uniprot-advanced-search-geobacillus-stearothermophilus.png)

Alternatively we could type this in the search bar:
- (protein_name:Laccase) AND (taxonomy_name:"Geobacillus stearothermophilus")

4. After hitting search we should see something like this:

![](figures/uniprot-advanced-search-geobacillus-stearothermophilus-results.png)

5. By clicking on the Download button and choosing 'FASTA' as the format we retrieve the sequences in a FASTA file. The file can be found under *data/* named *uniparc_protein_name_Laccase_AND_taxo_2026_02_16.fasta*.

6. Then we can proceed to visualization with Phylogenetic tree or Logo plots



