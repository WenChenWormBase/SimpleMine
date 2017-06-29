#!/usr/bin/perl -w

use strict;
use Ace;

print "This script create Gene-Disease-HumanOrtholog table based on disease association file and ortholog association file from the current WS release.\n";

my $acedbpath='/home/citace/WS/acedb/';
my $tace='/usr/local/bin/tace';

print "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) || die print "Connection failure: ", Ace->error;

#------------Build disease term and id hash-------

my %DO;  
my ($doTerm, $doName);

my $query="QUERY FIND DO_term";
my @doTermList = $db->find($query);
foreach $doTerm (@doTermList) {
     $doName = $doTerm->Name;
     $DO{$doTerm} = $doName;
     #print "$doTerm - $doName\n";
}
print scalar @doTermList, " disease ontology term found in WS.\n";

#--------------Done ---------------------------------------------------

#---------------- Get worm-human ortholog info -----------------
#open (ORTHO, "c_elegans.PRJNA13758.WS259.orthologs.txt")  || die "cannot open c_elegans.PRJNA13758.WSXXX.orthologs.txt!\n";
open (ORTHO, "/home/wen/simpleMine/ace_files/ce_orthologs.txt")  || die "cannot open c_elegans.PRJNA13758.WSXXX.orthologs.txt!\n";

my ($line, $spe);
my %orthoInfo;
my $geneid = "";
my $gName = "";
my $humanOrtho = "";
my $humanOrthoSingle;
my $humanid = "";
my $humanName = "";
my $evidence = "";
my $totalOrtho = 0;

while ($line=<ORTHO>) {
    chomp($line);
    if ($line eq "=") {
	if ($humanOrtho ne "") {
	    $orthoInfo{$geneid} = $humanOrtho;
	    $totalOrtho++;
	}
	$geneid = "";
	$gName = "";
	$humanOrtho = "";

    } elsif ($line =~ /^WBGene/) { #this line contains gene ID
	($geneid, $gName) = split /\s+/, $line;
    } elsif ($line =~ /^Homo sapiens/) { #this line contains 
	($spe, $humanid, $humanName, $evidence) = split /\t/, $line;
	
	if ($humanOrtho eq "") {
	    $humanOrtho = join '|', $humanid, $evidence;
	    
	} else {
	    $humanOrthoSingle = join '|', $humanid, $evidence;
	    $humanOrtho = join ',', $humanOrtho, $humanOrthoSingle;
	}
    }
}
close (ORTHO);
print "$totalOrtho C.elegans genes have human orthologs\n";

#--------------- Get Worm Disease info ---------------
#open (DAF, "disease_association.WS259.txt")  || die "cannot open disease_association.WSXXX.txt!\n";
open (DAF, "/home/wen/simpleMine/ace_files/disease_association.wb")  || die "cannot open disease_association.WSXXX.txt!\n";
my @tmp;
my %diseaseInfo;
my ($omim, $doEntry, $diseaseInformation);
my $totalDiseaseGene = 0;

while ($line=<DAF>) {
    chomp($line);
    next unless $line =~ /^WB/;
    @tmp = split /\s+/, $line;
    $geneid = $tmp[1];
    $gName = $tmp[2];
    $doTerm = $tmp[3];
    if ($tmp[5] eq "IEA") {
	$evidence = "Based on Sequence";
    } elsif ($tmp[5] eq "IMP") {
	$evidence = "Based on Experiment";
    } else {
	$evidence = "ERROR-No-DiseaseEvidence";
    }
    $omim = $tmp[6];
    if ($DO{$doTerm}) {
	$doName = $DO{$doTerm};
    } else {
	$doName = "N.A.";
    }
    
    $doEntry = join '|', $doName, $evidence, $doTerm, $omim;
    if ($diseaseInfo{$geneid}) {
	$diseaseInfo{$geneid} = join ',', $diseaseInfo{$geneid}, $doEntry;
    } else {
	$diseaseInfo{$geneid} = $doEntry;
	$totalDiseaseGene++;
    } 
}

close (DAF);
print "$totalDiseaseGene C.elegans genes have disease association.\n";

#-----------------Build Gene-Disease-HumanOrtholog table--------------
open (NAME, "/home/wen/simpleMine/ace_files/WBGeneIdentity.ace")  || die "cannot open /home/wen/simpleMine/ace_files/WBGeneIdentity.ace!\n";
    
open (OUT, ">GeneDiseaseHumanOrtholog.csv") || die "cannot open $!\n";
print OUT "WormBase Gene ID\tDisease Info\tHuman Ortholog\n";

while ($line = <NAME>) {

    #--- select C.elegans Gene only ----
    chomp($line);
    if ($line =~ /^Gene : /) {
	@tmp = split '"', $line;
	$geneid = $tmp[1];
    } elsif ($line =~ /^Public_name	 /) {
	@tmp = split '"', $line;
	$gName = $tmp[1];	
    } elsif ($line =~ /^Species	 /) {
	@tmp = split '"', $line;
	$spe = $tmp[1];

	if ($spe eq "Caenorhabditis elegans") {

	       if ($orthoInfo{$geneid}) {
		   $humanOrtho = ($orthoInfo{$geneid});
	       } else {
		   $humanOrtho = "N.A.";
	       }

	       if ($diseaseInfo{$geneid}) {
		   $diseaseInformation = $diseaseInfo{$geneid};
	       } else {
		   $diseaseInformation = "N.A.";
	       }

	       print OUT "$geneid\t$diseaseInformation\t$humanOrtho\n";
	    
	}
	
    }
    

}

close (OUT);
print "Done printing Gene-Disease-HumanOrtholog table.\n";
#-----------------------------------------------------------------------
$db->close();

