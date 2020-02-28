#!/usr/bin/perl

#-------------------Type out the purpose of the script-------------
print "This program create Gene Name table.\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneIdentity.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneTranscript.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneOperon.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WormPepLive.ace\n";
print "Output file: WBGeneName.csv\n\n";

print "Parse Gene names ...";

my ($line, $g, $cds, $pub_name, $merged_into, $status, $seq_name, $wormpep, $uniprot, $uniprotRef, $treefam, $refSeqRNA, $refSeqProtein, $other_name, $mol_name, $ope, $motif);
my @ceGenes;
my $ceID = 0;
my @tmp;
my %otherName;
my %pubName;
#my %cdsGene;
my %geneCDS;
my %wpExist;
my %wpMotif;
my %geneOPE;
my %geneUniprot;
my $geneUniPair;
my @otherNameList;
my ($w0, $w1); #for wormpep

#------------- Get transcript names ---------------------------------------------

open (CDS, "/home/wen/simpleMine/ace_files/WBGeneTranscript.ace") || die "can't open WBGeneTranscript.ace!";
while ($line=<CDS>) {
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
    if ($line =~ /^Protein/) {
	@tmp = split '"', $line;
	$wormpep = $tmp[1];
    } elsif ($line =~ /^Motif_homol/) {
	@tmp = split '"', $line;
	$motif = $tmp[1];
	if ($wpMotif{$wormpep}){
	    $wpMotif{$wormpep} = join ", ", $wpMotif{$wormpep}, $motif;
	    $wpExist{$wormpep}++;
	} else {
	    $wpMotif{$wormpep} = $motif;
	    $wpExist{$wormpep} = 1;
        }
    }
}
close (PEP);
#--------------------------- done getting WormPep names -------------------------------


#------------- Get operon names ---------------------------------
open (OPE, "/home/wen/simpleMine/ace_files/WBGeneOperon.ace") || die "can't open
 WBGeneOperon.ace!";
while ($line=<OPE>) {
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
while ($line =<IN>) {
    chomp($line);
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	$pub_name = "N.A.";
    } elsif ($line =~ /^Public_name/) {
	@tmp = split '"', $line;
	$pub_name = $tmp[1];
	$pubName{$pub_name} = $g;
    }
}
close (IN);
#------------ done  getting public name ------------------------------

open (IN, "/home/wen/simpleMine/ace_files/WBGeneIdentity.ace") || die "can't open WBGeneIdentity.ace!";
open (OUT, ">WBGeneName.csv") || die "cannot open WBGeneName.csv!\n";
open (CEG, ">AllCelegansGenes.txt") || die "cannot open AllCelegansGenes.txt!\n";

print OUT "WormBase Gene ID\tPublic Name\tWormBase Status\tSequence Name\tOther Name\tTranscript\tOperon\tWormPep\tProtein Domain\tUniprot\tReference Uniprot ID\tTreeFam\tRefSeq_mRNA\tRefSeq_protein\n";

$id = 0;
$id_a = 0;
$id_c = 0;
$oth = 0;

$line = <IN>;
while ($line =<IN>) {
    chomp($line);
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
	$merged_into = "N.A.";
	$other_name = "N.A.";

	$id++;
	print GENES "$id\t$g\n";
	$idGene{$g} = $id;
    } elsif ($line =~ /^Live/) {
	$status = "Live";
    } elsif ($line =~ /^Dead/) {	
	if ($merged_into ne "N.A.") {
	    $status = "Dead, merged into $merged_into";
	} else {
	    $status = "Dead";
	}

    } elsif ($line =~ /Caenorhabditis elegans/) {
	$ceGenes[$ceID] = $g;
	$ceID++;	

    } elsif ($line =~ /^Merged_into/) {
	@tmp = split '"', $line;
	$merged_into = $tmp[1];
    } elsif ($line =~ /^Sequence_name/) {
	@tmp = split '"', $line;
	$seq_name = $tmp[1];

	$id_a++;
	#print ALIAS "$id_a\t$seq_name\t$g\n";

    } elsif ($line =~ /^Molecular_name/) {
	@tmp = split '"', $line;
	$mol_name = $tmp[1];

	#get WormPep names		
	if ($wpExist{$mol_name}) {
	    $wormpep = $mol_name; 
	    if ($wpMotif{$mol_name}) {
	      $motif = join " - ", $wormpep, $wpMotif{$mol_name};
	    }	    
	    if ($allWormpep eq "N.A.") { #first protein product
		$allWormpep = $wormpep;
		$allMotif = $motif;
	    } else { #this gene has multiple protein product
		$allWormpep = join " \| ", $allWormpep, $wormpep;
		$allMotif = join " \| ", $allMotif, $motif;
	    }
	}

	$id_a++;
	#print ALIAS "$id_a\t$mol_name\t$g\n";
    } elsif ($line =~ /^Public_name/) {
	@tmp = split '"', $line;
	$pub_name = $tmp[1];
    } elsif ($line =~ /^Other_name/) {
	@tmp = split '"', $line;
	if ($publicName{$tm[1]}) {
	    #do not take it since this name is a public name for another gene
	} else {
	    if ($other_name eq "N.A.") {
		$other_name = $tmp[1];
	    } else {
		$other_name = join ", ",  $other_name, $tmp[1];
	    }
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
	
	print OUT "$g\t$pub_name\t$status\t$seq_name\t$other_name\t$cds\t$ope\t$allWormpep\t$allMotif\t$uniprot\t$uniprotRef\t$treefam\t$refSeqRNA\t$refSeqProtein\n";
    }
}

foreach $g (@ceGenes) {
    print CEG "$g\n";
}
print "$ceID C. elegans genes found.\n";
close (CEG);

close (IN);
#close (ALIAS);
#close (GENES);
#close (COMMON);
close (OUT);
print "Done.\n";


