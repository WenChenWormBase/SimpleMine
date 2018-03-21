#!/usr/bin/perl

#-------------------Type out the purpose of the script-------------
print "This program create Gene Name table.\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneIdentity.ace\n";
print "Input file: /home/wen/simpleMine/ace_files/WBGeneTranscript.ace\n";
print "Output file: WBGeneName.csv\n\n";

print "Parse Gene names ...";

my ($line, $g, $cds, $pub_name, $merged_into, $status, $seq_name, $wormpep, $uniprot, $treefam, $refSeqRNA, $refSeqProtein, $other_name, $mol_name);
my @ceGenes;
my $ceID = 0;
my @tmp;
my %otherName;
my %pubName;
#my %cdsGene;
my %geneCDS;
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

print OUT "WormBase Gene ID\tPublic Name\tWormBase Status\tSequence Name\tOther Name\tTranscript\tWormPep\tUniprot\tTreeFam\tRefSeq_mRNA\tRefSeq_protein\n";

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
	$uniprot = "N.A.";
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

	if ($mol_name =~ /^WP/) {
	    ($w0, $w1) = split ":", $mol_name;
		
	    if ($wormpep eq "N.A.") {
		$wormpep = $w1;
	    } else {
		$wormpep = join ", ",  $wormpep, $w1;
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
    }  elsif ($line =~ /UniProt/) {
	@tmp = split '"', $line;

	$id_a++;
	#print ALIAS "$id_a\t$tmp[5]\t$g\n";

	if ($uniprot eq "N.A.") {
	    $uniprot = $tmp[5];
	} else {
	    $uniprot = join ", ",  $uniprot, $tmp[5];
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

	print OUT "$g\t$pub_name\t$status\t$seq_name\t$other_name\t$cds\t$wormpep\t$uniprot\t$treefam\t$refSeqRNA\t$refSeqProtein\n";
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


