Simple Mine building pipeline

1. dumpWS_Gene_ace.sh
Dump out ace files for WBGeneName.csv

2. make_WBGene_table.pl
Create WBGeneName.csv and GeneHistory.csv from the ace files

3. makeRNAiAllelePheno.pl
Create RNAi and Allele phenotype table

4. makeGeneTissueLifeStage.pl
Create tissue and life stage table 

5. makeHomologDiseaseAssociation.pl
Create disease association and human homolog table

6. getGeneReference.pl
Create a table with WBPaper ID and PubMed ID for each gene

7. makeConfirmedInteraction.pl
Create a table for genes that showed Physical/Genetic/Regulatory interactions. There is an option for users to exlcude high-throughput interactions.

8. makeGeneAlleleTable.pl
Create a table for all sequenced alleles except those mapped to introns or are silent

9. makeGeneMap.pl
Create a table showing genetic map position, chromosome coordinates, and transcription factors

10. getGoTerm.pl
Create a table for Gene Ontology terms related to each gene

11. makeGeneDescription.pl
Create a table for Concise Description, Automatic Description and Expression Cluster Summary.

12. category_headers
A header file that divides search fields into categories. Need to update this file whenever search fields are added, deleted or renamed.

13. buildMultiSpecSimpleMine.sh
Create folders for species-specific SimpleMine

14. splitSpecies.pl
Split source data into species-specific folders built by buildMultiSpecSimpleMine.sh

15. buildSimpleMine.sh
A shell script that operate all above perl scripts

16. makeGeneNameSanitizer.pl
Create GeneNameHistory.csv file based on GeneHistory.csv and WBGeneName.ace. This is for the Gene Name Sanitizer tool. 
