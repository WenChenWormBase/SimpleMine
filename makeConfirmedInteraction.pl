#!/usr/bin/perl -w

use strict;

print "This script create Gene-ConfirmedInteraction table based on c_elegans.PRJNA13758.WSXXX.interactions.txt file from the current WS release. Interactions are confirmed when Genetic_interaction, Protein-Protein interaction and Predicted interactions all exist.\n";

my ($line, $g, $g1, $g2, $tmp_length, $genePair);
my @tmp;
my @intPreList; #for predicted interacting genes
my $i = 0; #this is the id for the interaction pairs that we will screen later.
my $j = 0; # this the id for the whole gene list to be printed
my @gPair;
my @sortedGP;
my %geneInt; #a list of interacting genes for each gene.
my %existGene; #this record if a gene is already in @geneList
my @geneList; #a list of genes that have interaction data
my %Predicted;
my %publicName;

#---- get public names for all genes ------------
open (IN, "/home/wen/simpleMine/ace_files/WBGeneIdentity.ace") || die "can't open WBGeneIdentity.ace!";
while ($line =<IN>) {
    chomp($line);
    if ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
    } elsif ($line =~ /^Public_name/) {
	@tmp = split '"', $line;
	$publicName{$g} = $tmp[1];
    }
}
close (IN);


#------ go through predicted interaction ------------
open (PRE, "/home/wen/simpleMine/ace_files/WSInt.ace")  || die "cannot open WSInt.ace!\n";
while ($line=<PRE>) {
    chomp($line);
    next unless ($line ne "");
    @tmp = ();
    
    @tmp = split '"', $line;
    if ($line =~ /^Interaction/) {
	#start a new pair
	$g1 = "",
	$g2 = "",
        @gPair = ();
	@sortedGP = ();
    } elsif ($line =~ /^Interactor_overlapping_gene/) {
	if ($g1 eq "") {
	    $g1 = $tmp[1];
	} else {
	    $g2 = $tmp[1];
	    @gPair = ($g1, $g2);
	    @sortedGP = sort @gPair;
	    $genePair = join '_', $sortedGP[0], $sortedGP[1];
	    if ($Predicted{$genePair}) {
		#do nothing, this pair is already recorded.
	    } else {
		$Predicted{$genePair} = 1;
		$intPreList[$i] = $genePair;
		$i++;
	    }
	}
    }
}
close (PRE);
print "$i pairs of genes found with Predicted interaction.\n";


open (INT, "/home/wen/simpleMine/ace_files/ce_interaction.txt")  || die "cannot open ce_interaction.txt!\n";
$i = 0;
my @intList;
my %Genetic; #Genetic interaction
my %PP; #Physical interaction
my %Reg; #Regulatory interaction
my %anyINT;#This record any gene pair interaction

while ($line=<INT>) {
    chomp($line);
    @tmp = ();
    @gPair = ();
    @sortedGP = ();
    @tmp = split /\t/, $line;
    $tmp_length = @tmp;
    next unless ($tmp_length == 11);
    
    #if ($tmp[2] eq "ProteinProtein") { # Protein-Protein interaction
    if ($tmp[1] eq "Physical") { # Protein-Protein interaction
	@gPair = ($tmp[5], $tmp[8]);
	@sortedGP = sort @gPair;
	$genePair = join '_', $sortedGP[0], $sortedGP[1];
	if ($PP{$genePair}) {
	    #do nothing, this pair is already recorded.
	} else {
	    $PP{$genePair} = 1;
	    if ($anyINT{$genePair}) {
		#do nothing, this pair is already recorded for the general list.
	    } else{
		$anyINT{$genePair} = 1;
		$intList[$i] = $genePair;
		$i++;
	    }
	}	
    #} elsif ($tmp[2] eq "Genetic_interaction") { # Genectic inertaction
    } elsif ($tmp[1] eq "Genetic") { # Genectic inertaction	
    	@gPair = ($tmp[5], $tmp[8]);
	@sortedGP = sort @gPair;
	$genePair = join '_', $sortedGP[0], $sortedGP[1];
	if ($Genetic{$genePair}) {
	    #do nothing, this pair is already recorded.
	} else {
	    $Genetic{$genePair} = 1;
	    if ($anyINT{$genePair}) {
		#do nothing, this pair is already recorded for the general list.
	    } else{
		$anyINT{$genePair} = 1;
		$intList[$i] = $genePair;
		$i++;
	    }
	}
    } elsif ($tmp[1] eq "Regulatory") { # Regulatory inertaction	
    	@gPair = ($tmp[5], $tmp[8]);
	@sortedGP = sort @gPair;
	$genePair = join '_', $sortedGP[0], $sortedGP[1];
	if ($Reg{$genePair}) {
	    #do nothing, this pair is already recorded.
	} else {
	    $Reg{$genePair} = 1;
	    if ($anyINT{$genePair}) {
		#do nothing, this pair is already recorded for the general list.
	    } else{
		$anyINT{$genePair} = 1;
		$intList[$i] = $genePair;
		$i++;
	    }	    
	}
    }    
}
close (INT);
#print "$i pairs of genes found with Protein-Protein interaction.\n";

my $evi;
my ($des1, $des2);
foreach $genePair (@intList) {
    #next unless (($Genetic{$genePair})||($PP{$genePair}));
    #this pair has confirmed interaction

    #record each gene into @geneList
    @tmp = ();
    @tmp = split '_', $genePair;
    foreach $g (@tmp) {
	if ($existGene{$g}) {
	    #skip, already recorded
	} else {
	    $geneList[$j] = $g;
	    $j++;
	    $existGene{$g} = 1;
	}	
    }

    #get interacting genes into %geneInt 
    $g1 = $tmp[0];
    $g2 = $tmp[1];
    next unless ($publicName{$g1});
    next unless ($publicName{$g2});

    #get types of interactions
    $evi = "";
    if ($PP{$genePair}) {
	$evi = "Physical";
    }
    if ($Genetic{$genePair}) {
	if ($evi eq "") {
	    $evi = "Genetic";
	} else {
	    #$evi = join " | ", $evi, "Genetic";
	    $evi = join ";", $evi, "Genetic";
	}
    }    
    if ($Reg{$genePair}) {
	if ($evi eq "") {
	    $evi = "Regulatory";
	} else {
	    #$evi = join " | ", $evi, "Regulatory";
	    $evi = join ";", $evi, "Regulatory";
	}
    }
    #done getting interaction types for each pair

    #$des1 = join "", $publicName{$g2}, "(", $evi, ")";;
    #$des2 = join "", $publicName{$g1}, "(", $evi, ")";;
    $des1 = join "\|", $publicName{$g2}, $evi;
    $des2 = join "\|", $publicName{$g1}, $evi;
   
    if ($geneInt{$g1}) {
	$geneInt{$g1} = join ", ", $geneInt{$g1}, $des1;
    } else {
	$geneInt{$g1} = $des1;
    }
    
    if ($geneInt{$g2}) {
	$geneInt{$g2} = join ", ", $geneInt{$g2}, $des2;
    } else {
	$geneInt{$g2} = $des2;
    }
    
}
print "$j genes found having confirmed interactions.\n";


my @Int3;
my @Int2;
my @Int1;
my ($allInt, $oneInt, $numEvi, $i1, $i2, $i3);
my @eviList;
my @sortedList;

open (OUT, ">ConfirmedInteraction.csv") || die "cannot open $!\n";
#print OUT "WormBase Gene ID\tConfirmed Interaction\n";
print OUT "WormBase Gene ID\tInteracting Gene\n";
foreach $g (@geneList) {
    if ($geneInt{$g}) {
	@tmp = ();
	@Int3 = ();
	@Int2 = ();
	@Int1 = ();
	$i1 = 0;
	$i2 = 0;
	$i3 =0;
	
	@tmp = split ", ", $geneInt{$g};
	@sortedList = ();
	@sortedList = sort { lc($a) cmp lc($b) } @tmp;
	foreach $oneInt (@sortedList) {
	    @eviList = ();
	    @eviList = split ";", $oneInt;
	    $numEvi = @eviList;
	    if ($numEvi == 1) {
		$Int1[$i1] = $oneInt;
		$i1++;
	    } elsif  ($numEvi == 2) {
		$Int2[$i2] = $oneInt;
		$i2++;
	    } elsif  ($numEvi == 3) {
		$Int3[$i3] = $oneInt;
		$i3++;
	    }
	    $allInt = join ", ", @Int3, @Int2, @Int1;
	}
	#print OUT "$g\t$geneInt{$g}\n";
	print OUT "$g\t$allInt\n";
    } else {
	print "ERROR: $g has no interaction!\n"
    }
}
close (OUT);
print "Done printing Gene-ConfirmedInteraction table.\n";
