#!/usr/bin/perl

use strict;


#-------------------Type out the purpose of the script-------------
print "This program creates Gene Name table.\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneIdentity.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneRemark.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneTranscript.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneOperon.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WormPepLive.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneSpe.ace\n";
print "Output file: WBGeneName.csv, GeneHistory.csv\n\n";

print "Parse Gene names ...";

my ($line, $g, $spe, $cds, $pub_name, $seq_name, $wormpep, $uniprot, $uniprotRef, $treefam, $refSeqRNA, $refSeqProtein, $other_name, $conflictName, $mol_name, $ope, $motif);
my @ceGenes;
my $ceID = 0;
my @tmp;
my %otherName;
my %pubName;
my %geneCDS;
my %geneSpe;
my %wpExist;
my %wpMotif;
my %geneOPE;
my %geneUniprot;
my $geneUniPair;
my @otherNameList;
my ($w0, $w1, $checkstring); #for wormpep
my %existPepMotif;

#------ Get Species information -------------------------
open (SPE, "/home/wen/simpleMine/ace_files/WBGeneSpe.ace") || die "can't open WBGeneSpe.ace!";
while ($line=<SPE>) {
    @tmp = ();
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
    } elsif ($line =~ /^Species/) {
	@tmp = split '"', $line;
	$spe = $tmp[1];
	$geneSpe{$g} = $spe;
    }
}
close (SPE);

#---- Done getting Species information -------------------



#------------- Get transcript names ---------------------------------------------

open (CDS, "/home/wen/simpleMine/ace_files/WBGeneTranscript.ace") || die "can't open WBGeneTranscript.ace!";
while ($line=<CDS>) {
    @tmp = ();  
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
    } elsif ($line =~ /^Corresponding_transcript/) {
	@tmp = split '"', $line;
	$cds = $tmp[1];
	#$cdsGene{$cds} = $g; 
	if ($geneCDS{$g}){
	    $geneCDS{$g} = join ", ", $geneCDS{$g}, $cds;
	} else {
	    $geneCDS{$g} = $cds;
        }
    }
}
close (CDS);
#--------------------------- done getting transcript names -------------------------------


#------------- Get WormPep names ---------------------------------------------
open (PEP, "/home/wen/simpleMine/ace_files/WormPepLive.ace") || die "can't open WormPepLive.ace!";
while ($line=<PEP>) {
    @tmp = ();      
    if ($line =~ /^Protein/) {
	@tmp = split '"', $line;
	$wormpep = $tmp[1];
    } elsif ($line =~ /^Motif_homol/) {
	@tmp = split '"', $line;
	$motif = $tmp[1];
	$checkstring = join "", $wormpep, $motif;
	next unless !($existPepMotif{$checkstring});
	if ($wpMotif{$wormpep}){
	    $wpMotif{$wormpep} = join "\|", $wpMotif{$wormpep}, $motif;
	    $wpExist{$wormpep}++;
	    $existPepMotif{$checkstring} = 1;
	} else {
	    $wpMotif{$wormpep} = $motif;
	    $wpExist{$wormpep} = 1;
	    $existPepMotif{$checkstring} = 1;
        }
    }
}
close (PEP);
#--------------------------- done getting WormPep names -------------------------------


#------------- Get operon names ---------------------------------
open (OPE, "/home/wen/simpleMine/ace_files/WBGeneOperon.ace") || die "can't open
 WBGeneOperon.ace!";
while ($line=<OPE>) {
    @tmp = ();  
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
    } elsif ($line =~ /^Contained_in_operon/) {
	@tmp = split '"', $line;
	$ope = $tmp[1];
	if ($geneOPE{$g}){
	    $geneOPE{$g} = join ", ", $geneOPE{$g}, $ope;
	} else {
	    $geneOPE{$g} = $ope;
        }
    }
}
close (OPE);
#--------------------------- done getting operon names -------------

#---------------- Get Public names ------------------------------------
open (IN, "/home/wen/simpleMine/ace_files/WBGeneIdentity.ace") || die "can't open WBGeneIdentity.ace!";
my %gPub; #locate public name when given a gene ID
 
while ($line =<IN>) {
    chomp($line);
    @tmp = ();  
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	$pub_name = "N.A.";
    } elsif ($line =~ /^Public_name/) {
	@tmp = split '"', $line;
	$pub_name = $tmp[1];
	$pubName{$pub_name} = $g;
	$gPub{$g} = $pub_name;
    }
}
close (IN);
#------------ done  getting public name ------------------------------

#---------------- Get Remark ------------------------------------
open (IN, "/home/wen/simpleMine/ace_files/WBGeneRemark.ace") || die "can't open WBGeneIdentity.ace!";
my $remark = "N.A.";
my %geneRemark;
while ($line =<IN>) {
    chomp($line);
    @tmp = ();  
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	$remark = "N.A.";
    } elsif ($line =~ /^Remark/) {
	@tmp = split '"', $line;
	if ($remark eq "N.A.") {
	    $remark = $tmp[1];
	} else {
	    $remark = join "...", $remark, $tmp[1];
	}
    } elsif ($line eq "") {
	if ($remark ne "N.A.") {
	    $geneRemark{$g} = $remark;
	}
    }
}
close (IN);
#------- done getting Remark -----

#------------ done  getting genes shareing the sequence name -------------------

open (IN, "/home/wen/simpleMine/ace_files/WBGeneIdentity.ace") || die "can't open WBGeneIdentity.ace!";
open (OUT, ">WBGeneName.csv") || die "cannot open WBGeneName.csv!\n";
open (CEG, ">AllCelegansGenes.txt") || die "cannot open AllCelegansGenes.txt!\n";
open (CON, ">GeneHistory.csv") || die "cannot open ConflictingNameGenes.txt!\n";

print CON "Gene Name\tPublic Name\tStatus\tSpecies\tIdentifier\tSplit Into\tMerged Into\tHistory\n";

print OUT "WormBase Gene ID\tPublic Name\tSpecies\tSequence Name\tOther Name\tTranscript\tOperon\tWormPep\tProtein Domain\tUniProt\tReference UniProt ID\tTreeFam\tRefSeq_mRNA\tRefSeq_protein\n";

my $id = 0;
my $id_a = 0;
my $id_c = 0;
my $oth = 0;
my ($allWormpep, $allMotif);
my %idGene;

#These are for gene history 
my @hisact;
my @deadGeneName;
my %deadGeneDes;
my $di = 0; #dead gene name list id
my ($status, $actTime, $actDate, $action, $history, $splitInto, $mergeInto);
my %nameStatus;
my %splitGene;
my %getAllNames; #Document all Public, Sequence and CGC names for each gene id;

$line = <IN>;
while ($line =<IN>) {
    chomp($line);
     @tmp = ();  
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	$pub_name = "N.A.";
	$seq_name = "N.A.";
	$wormpep  = "N.A.";
	$allWormpep = "N.A.";
	$motif = "N.A.";
	$allMotif = "N.A.";
	$uniprot = "N.A.";
	$uniprotRef = "N.A.";
	$treefam = "N.A.";
	$refSeqRNA = "N.A.";
	$refSeqProtein = "N.A.";
	$status = "N.A.";
	$other_name = "N.A.";
	$conflictName = "N.A.";
	$remark = "N.A.";
	$history = "N.A.";

	$splitInto = "";
	$mergeInto = "";
	
	$id++;
	print GENES "$id\t$g\n";
	$idGene{$g} = $id;

    } elsif ($line =~ /^Live/) {
	    $status = "Live";
    } elsif ($line =~ /^Dead/) {	
	    $status = "Dead";
    } elsif ($line =~ /^Suppressed/) {
	$status = "Suppressed";
    } elsif ($line =~ /^Split_into/) {	
        @tmp = split '"', $line;
	if ($splitGene{$g}) {
	    $splitGene{$g} = join ", ", $splitGene{$g}, "$gPub{$tmp[1]}\($tmp[1]\)";
	    $splitInto = join ", ", $splitInto, $tmp[1];
	} else {
	    $splitGene{$g} = "$gPub{$tmp[1]}\($tmp[1]\)";
	    $splitInto = $tmp[1];
	}
    } elsif ($line =~ /^Merged_into/) {	
        @tmp = split '"', $line;
	$mergeInto = $tmp[1];	
    } elsif ($line =~ /^Version_change/) {
	next unless !($line =~ /Imported/);
	@hisact = ();
	$actTime = "N.A.";
	$actDate = "N.A.";
	$action = "N.A.";
	@hisact = split /\s+/, $line;
	($actDate, $actTime) = split "_", $hisact[2];
	if ($hisact[4] eq "Other_name") {
	    $action = "Renamed on $actDate";
	} elsif ($hisact[4] eq "Merged_into") {
	    $action = "Merged into $hisact[5] on $actDate";
	} elsif ($hisact[4] eq "Split_into") {
	    $action = "Split into $hisact[5] on $actDate";
	} elsif (($hisact[4] eq "Killed")||($hisact[4] eq "Suppressed")||($hisact[4] eq "Resurrected")) {
	    $action = join " on ", $hisact[4], $actDate;
	}
	next unless ($action ne "N.A.");
	if ($history eq "N.A.") {
	    $history = $action;
	} else {
	    $history = join "; ", $history, $action;
	}
    } elsif ($line =~ /^Sequence_name/) {
	@tmp = split '"', $line;
	$seq_name = $tmp[1];
	$id_a++;
    } elsif ($line =~ /^Molecular_name/) {
	@tmp = split '"', $line;
	$mol_name = $tmp[1];

	#get WormPep names		
	if ($wpExist{$mol_name}) {
	    $wormpep = $mol_name; 
	    if ($wpMotif{$mol_name}) {
	      $motif = join "\|", $wormpep, $wpMotif{$mol_name};
	    }	    
	    if ($allWormpep eq "N.A.") { #first protein product
		$allWormpep = $wormpep;
		$allMotif = $motif;
	    } else { #this gene has multiple protein product
		$allWormpep = join ", ", $allWormpep, $wormpep;
		$allMotif = join ", ", $allMotif, $motif;
	    }
	}

	$id_a++;
	#print ALIAS "$id_a\t$mol_name\t$g\n";
    } elsif ($line =~ /^Public_name/) {
	@tmp = split '"', $line;
	$pub_name = $tmp[1];
    } elsif ($line =~ /^Other_name/) {
	@tmp = split '"', $line;
	if ($pubName{$tmp[1]}) {
	    #This name is a public name for another gene
	    $conflictName = $tmp[1];
	}
	if ($other_name eq "N.A.") {
	   $other_name = $tmp[1];
	} else {
	   $other_name = join ", ",  $other_name, $tmp[1];
	}
	
    }  elsif ($line =~ /UniProtAcc/) {

	@tmp = split '"', $line;
	
	if ($line =~ /UniProt_GCRP/) {
	    $uniprotRef = $tmp[5];
	}
	
	next unless (($line =~ /TrEMBL/)||($line =~ /SwissProt/));
	$geneUniPair = join "", $g, $tmp[5];

	if ($geneUniprot{$geneUniPair}) {
	    #do nothing
	} else {
	    $id_a++;
	    #print ALIAS "$id_a\t$tmp[5]\t$g\n";

	    if ($uniprot eq "N.A.") {
		$uniprot = $tmp[5];
	    } else {
		$uniprot = join ", ",  $uniprot, $tmp[5];
	    }
	    $geneUniprot{$geneUniPair} = 1;
	}
    }  elsif ($line =~ /TREEFAM/) {
	@tmp = split '"', $line;

	$id_a++;
	#print ALIAS "$id_a\t$tmp[5]\t$g\n";

	$treefam = $tmp[5];
    }  elsif (($line =~ /RefSeq/) && ($line =~ /mRNA/)) {
	@tmp = split '"', $line;

	$id_a++;
	#print ALIAS "$id_a\t$tmp[5]\t$g\n";
	
	if ($refSeqRNA eq "N.A.") {
	    $refSeqRNA = $tmp[5];
	} else {
	    $refSeqRNA = join ", ",  $refSeqRNA, $tmp[5];
	}
    }  elsif (($line =~ /RefSeq/) && ($line =~ /protein/)) {
	@tmp = split '"', $line;

	$id_a++;
	#print ALIAS "$id_a\t$tmp[5]\t$g\n";

	if ($refSeqProtein eq "N.A.") {
	    $refSeqProtein = $tmp[5];
	} else {
	    $refSeqProtein = join ", ",  $refSeqProtein, $tmp[5];
	}
    } elsif  ($line eq "") {
	if ($geneCDS{$g}) {
	    $cds = $geneCDS{$g};
	} else {
	    $cds = "N.A.";
	}

	if ($geneOPE{$g}) {
	    $ope = $geneOPE{$g};
	} else {
	    $ope = "N.A.";
	}

	if ($geneSpe{$g}) {
	    $spe = $geneSpe{$g};
	} else {
	    $spe = "N.A.";
	}

	if ($status eq "Live") {
	    print OUT "$g\t$pub_name\t$spe\t$seq_name\t$other_name\t$cds\t$ope\t$allWormpep\t$allMotif\t$uniprot\t$uniprotRef\t$treefam\t$refSeqRNA\t$refSeqProtein\n";

	    #next unless ($spe eq "Caenorhabditis elegans"); #only document C. elegans genes
	    #$ceGenes[$ceID] = $g;
	    #$ceID++;	

	    #if ($conflictName ne "N.A.") { #print ambiguous gene ids
		#$nameStatus{$conflictName} = "Ambiguous"; 
		#print CON "$conflictName\t\tPublic Name for $pubName{$conflictName}. May also refer to $pub_name\($g\). History of $g: $history\t\n";
	    #}

	    if ($splitGene{$g}) {
		print CON "$g\t$pub_name\t$status\t$spe\tLive - Split\t$splitInto\t$mergeInto\tSplit into $splitGene{$g}. History of $g: $history\t\n";		
		#print CON "$pub_name\t\tSplit into $splitGene{$g}. History of $pub_name\($g\): $history\t\n";
		#$nameStatus{$pub_name} = "Ambiguous";		
		#if (($seq_name ne $pub_name) && ($seq_name ne "N.A.")) {
		    #print CON "$seq_name\t\tSplit into $splitGene{$g}. History of $seq_name\($g\): $history\t\n";
		    #$nameStatus{$seq_name} = "Ambiguous";	
		#}
	    #} elsif ($shareSeqName{$seq_name}) {
		#two genes share the same sequence name
        	#print CON "$g\t$status\t$spe\tAmbiguous - Share Sequence\tTwo Genes share the same sequence: $shareSeqName{$seq_name}\t\n";
	        #print CON "$seq_name\t\tTwo Genes share the same sequence: $shareSeqName{$seq_name}.\t\n";		
		#$nameStatus{$seq_name} = "Ambiguous";
		#if ($pub_name ne "$seq_name") {
		    #print CON "$pub_name\t\tTwo Genes share the same sequence: $shareSeqName{$seq_name}.\t\n";
		    #$nameStatus{$pub_name} = "Ambiguous";
		#}		
	    } else {
		print CON "$g\t$pub_name\t$status\t$spe\tLive - Unique\t$splitInto\t$mergeInto\t$history\n"; #print valid gene ids
	    }

	    #---------------------------------------
	    #$getAllNames{$g} = join ", ", $pub_name, $seq_name, $other_name;
	    #-------------------------------------------
	    
	} else { #print obsolete gene ids
	    #if (($pub_name ne "N.A.")&& ($spe eq "Caenorhabditis elegans")) {
		#$deadGeneName[$di] = $pub_name;
		#$di++;
		#$deadGeneDes{$pub_name} = "$status gene. History of $pub_name\($g\): $history"
	    #}
	    print CON "$g\t$pub_name\t$status\t$spe\tObsolete - $status\t$splitInto\t$mergeInto\t$history";
	    if ($geneRemark{$g}) {
	       print CON ". Remark: $geneRemark{$g}\n";
	    } else  {
	       print CON "\n";
	    }
	}

    }
}
close (IN);
close (OUT);
#print "Found $di dead gene with public names.\n";

#my @allNames;
#my $genename;
#foreach $g (@ceGenes) {
    #print CEG "$g\n";
    #next unless ($getAllNames{$g});
    #@allNames = ();
    #$genename = "";
    #@allNames = split ", ", $getAllNames{$g};
    #foreach $genename (@allNames) {
	#next unless ($genename ne "N.A.");
	#if ($nameStatus{$genename}) {	    
	    #skip, this name was already printed 
	#} else {
	    #print CON "$genename\t$g\t\t\n";
	    #$nameStatus{$genename} = "Printed as Valid";
	#} 
    #}
#}
#print "$ceID C. elegans genes found.\n";

#foreach $genename (@deadGeneName) {
    #if ($nameStatus{$genename}) {
	#this name was already used by a Live gene
    #} else {
	#print CON "$genename\t\t\t$deadGeneDes{$genename}\n";
	#$nameStatus{$genename} = "Printed as Obsolete";
    #}
#}

close (CEG);
close (CON);
print "Done.\n";


#---------- print the rest of valide gene names -------
#open (IN, "/home/wen/simpleMine/ace_files/WBGeneName.ace") || die "can't open WBGeneName.ace!";
#my ($genename, $nameMatch);
#my %nameGene;
#my @validGeneName;
#my $vi = 0;
#while ($line =<IN>) {
#    chomp($line);
#    @tmp = ();  
#    if ($line =~ /^Gene_name/) {
#	@tmp = split '"', $line;
#	$genename = $tmp[1];
#	$g = "N.A.";
#	$nameMatch = 0;
#    } elsif ($line =~ /WBGene/) {
#	@tmp = split '"', $line;
#	if ($g eq "N.A.") {
#	    next unless ($nameStatus{$tmp[1]}); #only consider live gene ids
#	    $g = $tmp[1];
#	    $nameMatch++;
#	} elsif ($g eq $tmp[1])  {
#	    #skip, this gene id is already documented
#	} else {
#	    #this is a name with multiple match
#	    if ($nameStatus{$genename}) {
#		#This name was already documented as an ambigous name. 
#	    } else {
		#print "$genename has multiple id match but not recorded as an ambiguous name.\n";
#	    }
#	} 
#    } elsif ($line eq "") {
#	if ($nameMatch == 1) {
	    #print CON "$genename\t$g\t\t\n";
	    #$nameGene{$genename} = $g;
	    #$validGeneName{$vi} = $genename;
	    #$vi++;
#	}
#	$genename = "";
#    }
#}
#print "Printed $vi names that match only one WormBase gene id.\n";
#close (IN);


