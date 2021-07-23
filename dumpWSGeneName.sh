#!/bin/csh
cd /home/wen/simpleMine/
setenv ACEDB /home/citace/WS/acedb/
## from Wen
/usr/local/bin/tace -tsuser 'wen' <<END_TACE 
QUERY FIND Gene_name
show -a -f ace_files/WBGeneName.ace
quit
END_TACE

