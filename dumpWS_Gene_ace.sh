#!/bin/csh

#------- prepare disease association and homolog files -----------------
cp /home/citace/WS/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS*.orthologs.txt.gz /home/wen/simpleMine/ace_files/ce_orthologs.txt.gz
cd /home/wen/simpleMine/ace_files/
gunzip ce_orthologs.txt.gz
cp /home/citace/WS/ONTOLOGY/disease_association.WS*.wb /home/wen/simpleMine/ace_files/disease_association.wb

#-------------prepare ace files for PCL file generation---------------------
cd /home/wen/simpleMine/
setenv ACEDB /home/citace/WS/acedb/
## from Wen
/usr/local/bin/tace -tsuser 'wen' <<END_TACE
QUERY FIND Gene WBGene*
show -a -t Identity -f ace_files/WBGeneIdentity.ace
show -a -t Species -f ace_files/WBGeneSpe.ace
show -a -t Molecular_info -f ace_files/WBGeneTranscript.ace
quit
END_TACE
