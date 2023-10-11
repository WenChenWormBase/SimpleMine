#!/bin/csh

#------- prepare flat files -----------------
cd /home/wen/simpleMine/ace_files/
rm /home/wen/simpleMine/ace_files/*
#cp /home/wen/simpleMine/bin/OpenBiosystemsRNAiLibrary.csv /home/wen/simpleMine/ace_files/.
cp /home/citace/WS/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS*.annotations.gff3.gz /home/wen/simpleMine/ace_files/ce.gff3.gz
gunzip /home/wen/simpleMine/ace_files/ce.gff3.gz
cp /home/citace/WS/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS*.pcr_product2gene.txt.gz /home/wen/simpleMine/ace_files/ce_product2gene.txt.gz
gunzip /home/wen/simpleMine/ace_files/ce_product2gene.txt.gz
cp /home/citace/WS/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS*.orthologs.txt.gz /home/wen/simpleMine/ace_files/ce_orthologs.txt.gz
gunzip /home/wen/simpleMine/ace_files/ce_orthologs.txt.gz
cp /home/citace/WS/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS*.interactions.txt.gz /home/wen/simpleMine/ace_files/ce_interaction.txt.gz
gunzip /home/wen/simpleMine/ace_files/ce_interaction.txt.gz
cp /home/citace/WS/ONTOLOGY/gene_association.WS*.wb /home/wen/simpleMine/ace_files/gene_association.wb
cp /home/citace/WS/ONTOLOGY/gene_ontology.WS*.obo /home/wen/simpleMine/ace_files/gene_ontology.obo
#cp /home/citace/WS/ONTOLOGY/disease_association.WS*.wb /home/wen/simpleMine/ace_files/disease_association.wb
cat /home/wen/AutoDescription/ecSummary/*.csv > /home/wen/simpleMine/ace_files/allECreg.csv



#-------------prepare ace files for PCL file generation---------------------
cd /home/wen/simpleMine/
setenv ACEDB /home/citace/WS/acedb/
## from Wen
/usr/local/bin/tace -tsuser 'wen' <<END_TACE 
QUERY FIND Gene WBGene*; Remark = *
show -a -t Remark -f ace_files/WBGeneRemark.ace
QUERY FIND Gene WBGene*;
show -a -t Identity -f ace_files/WBGeneIdentity.ace
show -a -t Species -f ace_files/WBGeneSpe.ace
show -a -t Structured_description -f ace_files/WBGeneDescription.ace
show -a -t Molecular_info -f ace_files/WBGeneTranscript.ace
show -a -t Contained_in_operon -f ace_files/WBGeneOperon.ace
show -a -t Biotype -f ace_files/WBGeneBiotype.ace
QUERY FIND Clone PCR_product = *;
show -a -f ace_files/PCR_product.ace
QUERY FIND SO_term;
show -a -t SO_name -f ace_files/SO_name.ace
QUERY FIND Protein; WormPep AND Live 
show -a -t Motif_homol -f ace_files/WormPepLive.ace
QUERY FIND Variation Variation_type = *allele; Sequenced
show -a -f ace_files/SequencedAllele.ace
QUERY FIND Interaction High_throughput 
show -a -t High_throughput -f ace_files/HTPInt.ace
QUERY FIND Gene Map_info = *
show -a -t Map_info -f ace_files/WBGeneMap.ace
QUERY FIND Gene_name
show -a -f ace_files/WBGeneName.ace
quit
END_TACE
