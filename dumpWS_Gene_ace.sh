#!/bin/csh
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
