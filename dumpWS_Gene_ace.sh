#!/bin/csh

#------- prepare disease association and homolog files -----------------
cd /home/wen/simpleMine/ace_files/
rm /home/wen/simpleMine/ace_files/*
cp /home/citace/WS/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS*.orthologs.txt.gz /home/wen/simpleMine/ace_files/ce_orthologs.txt.gz
gunzip /home/wen/simpleMine/ace_files/ce_orthologs.txt.gz
cp /home/citace/WS/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS*.interactions.txt.gz /home/wen/simpleMine/ace_files/ce_interaction.txt.gz
gunzip /home/wen/simpleMine/ace_files/ce_interaction.txt.gz
#cp /home/citace/WS/ONTOLOGY/disease_association.WS*.wb /home/wen/simpleMine/ace_files/disease_association.wb
cat /home/wen/AutoDescription/ecSummary/*.csv > /home/wen/simpleMine/ace_files/allECreg.csv

#-------------prepare ace files for PCL file generation---------------------
cd /home/wen/simpleMine/
setenv ACEDB /home/citace/WS/acedb/
## from Wen
/usr/local/bin/tace -tsuser 'wen' <<END_TACE
QUERY FIND Gene WBGene*; Status = Live 
show -a -t Identity -f ace_files/WBGeneIdentity.ace
show -a -t Species -f ace_files/WBGeneSpe.ace
show -a -t Structured_description -f ace_files/WBGeneDescription.ace
show -a -t Molecular_info -f ace_files/WBGeneTranscript.ace
show -a -t Contained_in_operon -f ace_files/WBGeneOperon.ace
QUERY FIND Protein; WormPep AND Live 
show -a -t Motif_homol -f ace_files/WormPepLive.ace
QUERY FIND Variation Variation_type = *allele; Sequenced
show -a -f ace_files/SequencedAllele.ace
QUERY FIND Interaction Predicted AND Log_likelihood_score > 1.5 
show -a -t Interactor_overlapping_gene -f ace_files/WSInt.ace
QUERY FIND Gene Map_info = *
show -a -t Map_info -f ace_files/WBGeneMap.ace
quit
END_TACE
