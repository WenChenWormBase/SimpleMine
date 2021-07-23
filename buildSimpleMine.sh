cd /home/wen/simpleMine/
/home/wen/simpleMine/bin/buildMultiSpecSimpleMine.sh
mkdir /home/wen/simpleMine/sourceFile/
mkdir /home/wen/simpleMine/sourceFile/obsolete/
mkdir /home/wen/simpleMine/GeneNameSanitizer/

cd /home/wen/simpleMine/sourceFile/
/home/wen/simpleMine/bin/dumpWS_Gene_ace.sh
/home/wen/simpleMine/bin/make_WBGene_table.pl
/home/wen/simpleMine/bin/getGeneReference.pl
/home/wen/simpleMine/bin/getGoTerm.pl
/home/wen/simpleMine/bin/makeRNAiAllelePheno.pl
/home/wen/simpleMine/bin/makeGeneTissueLifeStage.pl
/home/wen/simpleMine/bin/makeHomologDiseaseAssociation.pl
/home/wen/simpleMine/bin/makeConfirmedInteraction.pl
/home/wen/simpleMine/bin/makeGeneAlleleTable.pl
/home/wen/simpleMine/bin/makeGeneMap.pl
/home/wen/simpleMine/bin/makeGeneDescription.pl
mv RNAiAllelePheno.csv obsolete/.
cp WBGeneName.csv 1_WBGeneName.csv
mv GeneticMapPosition.csv 2_GeneticMapPosition.csv
mv RNAiAllelePhenoObserved.csv 3_RNAiAllelePhenoObserved.csv
mv GeneAllele.csv 4_GeneAllele.csv
mv ConfirmedInteraction.csv 5_ConfirmedInteraction.csv
mv GeneTissueLifeStage.csv 6_GeneTissueLifeStage.csv
mv GeneDiseaseHumanOrtholog.csv 7_GeneDiseaseHumanOrtholog.csv
mv GeneReference.csv 8_GeneReference.csv
mv GeneOntologyAssociation.csv 9_GeneOntologyAssociation.csv
mv GeneDescription.csv 10_GeneDescription.csv

/home/wen/simpleMine/bin/makeGeneNameSanitizer.pl
mv GeneNameHistory.csv /home/wen/simpleMine/GeneNameSanitizer/.

/home/wen/simpleMine/bin/splitSpecies.pl
