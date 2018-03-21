Simple Mine building pipeline

1. dumpWS_Gene_ace.sh
dump out ace files for WBGeneName.csv

2. make_WBGene_table.pl
create WBGeneName.csv from the ace files

3. makeRNAiAllelePheno.pl
create RNAi and Allele phenotype table

4. makeGeneTissueLifeStage.pl
create tissue and life stage table 

5. makeHomologDiseaseAssociation.pl
create disease association and human homolog table

6. getGeneReference.pl
Create a table with WBPaper ID and PubMed ID for each gene

7. makeConfirmedInteraction.pl
Create a table for genes that showed Physical/Genetic/Regulatory interactions

8. makeGeneAlleleTable.pl
Create a table for all sequenced alleles except those mapped to introns or are silent

9. makeGeneMap.pl
Create a table showing chromosome and chromosomal position

10. buildSimpleMine.sh
A shell script that operate all above perl scripts
