#!/usr/bin/perl -w

use strict;

print "This script create Gene-ConfirmedInteraction table based on c_elegans.PRJNA13758.WSXXX.interactions.txt file from the current WS release. Interactions are confirmed when Genetic_interaction, Protein-Protein interaction and Predicted interactions all exist.\n";

my ($line, $g, $g1, $g2, $tmp_length, $genePair, $i);
my @tmp;
my %HTP; #high-throughput interaction objects
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
open (HTP, "/home/wen/simpleMine/ace_files/HTPInt.ace")  || die "cannot open l_HTP_Int.ace!\n";
$line=<HTP>;
$i = 0;
while ($line=<HTP>) {
    next unless ($line =~ /^Interaction/);
    @tmp = split '"', $line;
    $HTP{$tmp[1]} = 1;
    $i++;
}
close (HTP);
print "$i high-throughput interaction objects found.\n";



#----- Process this twice, the first time for all int, the second time for LTP only ----
my $c = 0; #this is for the two round of reading interaction file 
my $j = 0; # this the id for the whole gene list to be printed
my @gPair;
my @sortedGP;
my %geneInt; #a list of interacting genes for each gene.
my %existGene; #this record if a gene is already in @geneList
my @geneList; #a list of genes that have interaction data
my @intList;
my %Genetic; #Genetic interaction
my %PP; #Physical interaction
my %Reg; #Regulatory interaction
my %anyINT;#This record any gene pair interaction
my $htp;
my $evi;
my ($des1, $des2);
my @Int3;
my @Int2;
my @Int1;
my ($allInt, $oneInt, $numEvi, $i1, $i2, $i3);
my @eviList;
my @sortedList;

my %gAllInt; #all interactions for each gene
my %gLTPInt; #low-throughput interactions for each gene

while ($c < 2) {

    %geneInt = (); #a list of interacting genes for each gene.
    #%existGene = (); #this record if a gene is already in @geneList
    #@geneList = ();     
    @intList = ();
    %Genetic = (); #Genetic interaction
    %PP = (); #Physical interaction
    %Reg = (); #Regulatory interaction
    %anyINT = ();#This record any gene pair interaction
    
    open (INT, "/home/wen/simpleMine/ace_files/ce_interaction.txt")  || die "cannot open ce_interaction.txt!\n";
    $i = 0;
    while ($line=<INT>) {
	chomp($line);
	@tmp = ();
	@gPair = ();
	@sortedGP = ();
	@tmp = split /\t/, $line;
	$tmp_length = @tmp;
	next unless ($tmp_length == 11);

	$htp = 0;
	if ($c == 1) { #the second round to exclude HTP interactions 
	    if ($HTP{$tmp[0]}) {
		#this is a high-throughput object
		$htp = 1;
	    }
	}
	next unless ($htp == 0);
    
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

    foreach $genePair (@intList) {

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
	    if ($c == 0) {
		$gAllInt{$g} = $allInt;
	    } else {
		$gLTPInt{$g} = $allInt;
	    }
	} else {
	    #print "ERROR: $g has no interaction!\n"
	}
    }
    $c++;
}
print "$j genes found having confirmed interactions including high-throughput.\n";


my ($ATP, $LTP);
open (OUT, ">ConfirmedInteraction.csv") || die "cannot open $!\n";
print OUT "WormBase Gene ID\tInteracting Gene\tInteracting Gene Exclude High-throughput\n";
foreach $g (@geneList) {
    $ATP = "N.A.";
    $LTP = "N.A.";
    if ($gAllInt{$g}) {
	$ATP = $gAllInt{$g};
    }
    if ($gLTPInt{$g}) {
	$LTP = $gLTPInt{$g};
    }    
    print OUT "$g\t$ATP\t$LTP\n";
}
close (OUT);
print "Done printing Gene-ConfirmedInteraction table.\n";

