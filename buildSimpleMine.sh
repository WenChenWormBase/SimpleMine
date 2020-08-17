mkdir /home/wen/simpleMine/sourceFile/
mkdir /home/wen/simpleMine/sourceFile/obsolete/
mkdir /home/wen/simpleMine/GeneNameSanitizer/ 
mkdir multiSpeSimpleMine/b_malayi
mkdir multiSpeSimpleMine/c_briggsae
mkdir multiSpeSimpleMine/c_brenneri
mkdir multiSpeSimpleMine/c_elegans
mkdir multiSpeSimpleMine/c_japonica
mkdir multiSpeSimpleMine/c_remanei
mkdir multiSpeSimpleMine/p_pacificus
mkdir multiSpeSimpleMine/o_volvulus
mkdir multiSpeSimpleMine/s_ratti
mkdir multiSpeSimpleMine/mix

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
mv GeneOntologyAssociation.csv 9_1_GeneOntologyAssociation.csv
mv GeneDescription.csv 9_2_GeneDescription.csv

/home/wen/simpleMine/bin/splitSpecies.pl
