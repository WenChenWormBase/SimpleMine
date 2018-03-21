mkdir /home/wen/simpleMine/sourceFile/
cd /home/wen/simpleMine/sourceFile/
mkdir obsolete/
/home/wen/simpleMine/bin/dumpWS_Gene_ace.sh
/home/wen/simpleMine/bin/make_WBGene_table.pl
/home/wen/simpleMine/bin/getGeneReference.pl
/home/wen/simpleMine/bin/makeRNAiAllelePheno.pl
/home/wen/simpleMine/bin/makeGeneTissueLifeStage.pl
/home/wen/simpleMine/bin/makeHomologDiseaseAssociation.pl
/home/wen/simpleMine/bin/makeConfirmedInteraction.pl
/home/wen/simpleMine/bin/makeGeneAlleleTable.pl
/home/wen/simpleMine/bin/makeGeneMap.pl

mv RNAiAllelePheno.csv obsolete/.
cp WBGeneName.csv 1_WBGeneName.csv
mv GeneticMapPosition.csv 2_GeneticMapPosition.csv
mv RNAiAllelePhenoObserved.csv 3_RNAiAllelePhenoObserved.csv
mv GeneAllele.csv 4_GeneAllele.csv
mv ConfirmedInteraction.csv 5_ConfirmedInteraction.csv
mv GeneTissueLifeStage.csv 6_GeneTissueLifeStage.csv
mv GeneDiseaseHumanOrtholog.csv 7_GeneDiseaseHumanOrtholog.csv
mv GeneReference.csv 8_GeneReference.csv

