#!/usr/bin/perl -w

use strict;

print "This script locate Genetic map positions, chromosome coordinates, and transcription factors table based on the current WS release via ace files and gff3.\n";

my ($line, $g, $chr, $p, $err, $map, $tmp_length);
my @tmp;
my @pos;
my $i = 0;
my $totalGenes;
#my @geneMapList; #all genes with genetic map positions
my @geneList; #all genes with genetic or chromosome map positions
my %geneMap; #genetic map positions
my %existGene; #genes that were already recorded by this script
open (IN, "/home/wen/simpleMine/ace_files/WBGeneMap.ace") || die "can't open WBGeneMap.ace!";

while ($line = <IN>) {
    chomp ($line);
    next unless ($line ne "");

    @tmp = ();
    @pos = ();
    if ($line=~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	next unless ($g =~ /^WBGene/);
	#$geneMapList[$i] = $g;
	$geneList[$i] = $g;
	$existGene{$g} = 1;
	$i++;
	$chr = "";
	$p = "";
	$err = "";
    } elsif ($line =~ /^Map/) {
	next unless ($line =~ /"/);
	@tmp = split '"', $line;
	$chr = $tmp[1];

	if ($line =~ /Position/) {
	    @pos = split /\s+/, $line;
	    $p = sprintf("%.2f", $pos[3]);
	    if ($line =~ /Error/) {
		$err =  sprintf("%.2f", $pos[5]);
	    } 
	}

	if ($p eq "") {
	    #print OUT "$g\t$chr\n";
	    $geneMap{$g} = "$chr";
	} elsif ($err eq "") {
	    #print OUT "$g\t$chr $p\n";
	    $geneMap{$g} = "$chr $p";
	} else {
	    $geneMap{$g} = "$chr $p +\/- $err";
	}
    } elsif ($line =~ /^Interpolated_map_position/) {
	@tmp = split '"', $line;
	$chr = $tmp[1];

	@tmp = ();
	@tmp = split /\s+/, $line;
	
	$p = $tmp[2];
	#print OUT "$g\t$chr $p\n";
        $geneMap{$g} = "$chr $p";	
    }
}    
$totalGenes = $i;
close (IN);
print "Found genetic map positions for $totalGenes genes.\n";


my ($strand, $start, $lateStart, $stop, $stuff_length);
my @stuff;
my $totalGffGenes = 0;
my %geneStart;
my %geneStop;
my %geneChr;
my %geneStrand;
my %geneLateStart; #the latest start location of all transcripts of the gene
my %existTF;
my @data;
my ($dataline, $col, $bStart, $bStop, $dis, $tf, $gTF, $tfList);
my @sortTF;

open (IN, "/home/wen/simpleMine/ace_files/ce.gff3") || die "can't find ce.gff3!";

open (IDATA, ">/home/wen/simpleMine/ace_files/I_TF_bindingSite.csv") || die "cannot open I_TF_bindingSite.csv!\n";
open (IIDATA, ">/home/wen/simpleMine/ace_files/II_TF_bindingSite.csv") || die "cannot open II_TF_bindingSite.csv!\n";
open (IIIDATA, ">/home/wen/simpleMine/ace_files/III_TF_bindingSite.csv") || die "cannot open III_TF_bindingSite.csv!\n";
open (IVDATA, ">/home/wen/simpleMine/ace_files/IV_TF_bindingSite.csv") || die "cannot open IV_TF_bindingSite.csv!\n";
open (VDATA, ">/home/wen/simpleMine/ace_files/V_TF_bindingSite.csv") || die "cannot open V_TF_bindingSite.csv!\n";
open (XDATA, ">/home/wen/simpleMine/ace_files/X_TF_bindingSite.csv") || die "cannot open X_TF_bindingSite.csv!\n";
open (MDATA, ">/home/wen/simpleMine/ace_files/Mt_TF_bindingSite.csv") || die "cannot open Mt_TF_bindingSite.csv!\n";

while ($line = <IN>) {
    chomp ($line);
    @tmp = ();
    @tmp = split /\t/, $line;
    $tmp_length = @tmp;
    next unless ($tmp_length > 8);  
    if (($tmp[1] eq "WormBase")&&($tmp[2] eq "gene")) {
  	$chr = $tmp[0];
    	$start = $tmp[3];
    	$stop = $tmp[4];
    	$strand = $tmp[6];
    	@stuff = ();
    	@stuff = split /[\:\;]/, $tmp[8];
    	$g = $stuff[1];
	if ($existGene{$g}) {
	    #do nothing, this gene is already recorded at genetic map position
	} else {
	    print "$g has chromosome coordinates but not genetic map position!\n";
	    $totalGenes++;
	    $geneList[$totalGenes] = $g; #add this gene to the total list
	}
    	$geneChr{$g} = $chr;
    	$geneStart{$g} = $start;
	$geneStop{$g} = $stop;
    	$geneStrand{$g} = $strand;    
    	$totalGffGenes++; #document the number of genes on GFF3

        #document the late start site for each gene  
	if ($geneLateStart{$g}) {
	    $lateStart = $geneLateStart{$g};
	    if ($strand eq "-") {
		if ($lateStart < $stop) {
		    $geneLateStart{$g} = $lateStart;
		}
	    } else {#gene is on postive strand
		if ($lateStart > $start) {
		    $geneLateStart{$g} = $lateStart;
		}
	    }
	} else {
	    if ($strand eq "-") {#gene is on negative strand
		$geneLateStart{$g} = $stop;		
	    } else {#gene is on postive strand
		$geneLateStart{$g} = $start;
	    }
	}	
    } elsif (($tmp[1] eq "WormBase")&&($tmp[2] eq "mRNA")) {	
    	@stuff = ();
    	@stuff = split /[\:\;]/, $tmp[8];
    	$g = $stuff[3];	
	$strand = $tmp[6];

	#document the late start site for each gene
	if ($strand eq "-") {#transcript is on negative strand
	    $lateStart = $tmp[4];
	    if ($geneLateStart{$g}) {
		next unless  ($lateStart < $geneLateStart{$g});
		$geneLateStart{$g} = $lateStart;
	    } else {
		$geneLateStart{$g} = $lateStart;
	    }
	} else {#gene is on postive strand
	    $lateStart = $tmp[3];
	    if ($geneLateStart{$g}) {
		next unless  ($lateStart > $geneLateStart{$g});
		$geneLateStart{$g} = $lateStart;
	    } else {
		$geneLateStart{$g} = $lateStart;
	    }
        }
	
    } elsif ($tmp[2] =~ /TF_binding_site/) {
	@stuff = (); 
	$tf = "";
	if ($tmp[8] =~ /tf_name/) {
	    @stuff = split '=', $tmp[8];
	    $stuff_length = @stuff;
	    next unless ($stuff_length > 3);
            $tf = $stuff[3];	
	} elsif ($tmp[8] =~ /_Tf_name/) {
	    @stuff = split /[\:\;]/, $tmp[8];
	    $stuff_length = @stuff;
	    next unless ($stuff_length > 3);
 	    $tf = $stuff[3];
	} elsif ($tmp[8] =~ /ChIP-Seq/) {
            @stuff = split /\s/, $tmp[8];
	    $stuff_length = @stuff;
	    next unless ($stuff_length > 5);
            $tf = $stuff[5];
	} else {
		#print "ERROR! Cannot locate TF for $line!\n";
	}
	if ($tf ne "") {
		$chr = $tmp[0];
        	$start = $tmp[3];
        	$stop = $tmp[4];
        	$strand = $tmp[6];
		if ($chr eq "I") {
		    print IDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} elsif ($chr eq "II") {
		    print IIDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} elsif ($chr eq "III") {
		    print IIIDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} elsif ($chr eq "IV") {
		    print IVDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} elsif ($chr eq "V") {
		    print VDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} elsif ($chr eq "X") {
		    print XDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} elsif ($chr eq "MtDNA") {
		    print MDATA "$chr\t$strand\t$start\t$stop\t$tf\n";
		} else {
		    print "Error! Cannot find chromosome for $tf in $line!\n";
		}		
		#print  DATA "$chr\t$strand\t$start\t$stop\t$tf\n";
	}
    }
}
close (IN);
close (IDATA);
close (IIDATA);
close (IIIDATA);
close (IVDATA);
close (VDATA);
close (XDATA);
close (MDATA);
print "Found $totalGffGenes genes on GFF3. There are $totalGenes with either genetic map position or chromosome coordinates.\n";

open (OUT, ">GeneticMapPosition.csv") || die "cannot open GeneticMapPosition.csv!\n";
print OUT "WormBase Gene ID\tGenetic Map Position\tChromosome Coordinates\tTranscription Factors\n";

foreach $g (@geneList) {
    if ($geneMap{$g}) {
	$map = $geneMap{$g};
    } else {
	$map = "N.A.";
    }

    if ($geneChr{$g}) { #this gene is on GFF3
	
	$chr = $geneChr{$g};
	$strand = $geneStrand{$g};
	$start = $geneStart{$g};
	$stop = $geneStop{$g};
	$lateStart = $geneLateStart{$g};
	
	if ($chr eq "I") {
	    open (DATA, "/home/wen/simpleMine/ace_files/I_TF_bindingSite.csv") || die "Can't open I_TF_bindingSite.csv!";
	} elsif  ($chr eq "II") {
	    open (DATA, "/home/wen/simpleMine/ace_files/II_TF_bindingSite.csv") || die "Can't open II_TF_bindingSite.csv!";
	} elsif  ($chr eq "III") {
	    open (DATA, "/home/wen/simpleMine/ace_files/III_TF_bindingSite.csv") || die "Can't open III_TF_bindingSite.csv!";
	} elsif  ($chr eq "IV") {
	    open (DATA, "/home/wen/simpleMine/ace_files/IV_TF_bindingSite.csv") || die "Can't open IV_TF_bindingSite.csv!";
	} elsif  ($chr eq "V") {
	    open (DATA, "/home/wen/simpleMine/ace_files/V_TF_bindingSite.csv") || die "Can't open V_TF_bindingSite.csv!";
	} elsif  ($chr eq "X") {
	    open (DATA, "/home/wen/simpleMine/ace_files/X_TF_bindingSite.csv") || die "Can't open X_TF_bindingSite.csv!";
	} elsif  ($chr eq "MtDNA") {
	    open (DATA, "/home/wen/simpleMine/ace_files/Mt_TF_bindingSite.csv") || die "Can't open Mt_TF_bindingSite.csv!";
	} else {
	    print "ERROR! Cannot find chromosome information for $g!\n";
	}

        $tfList = "";
        while ($dataline = <DATA>) {
                chomp ($dataline);
		my $foundTF = 0;
                @data = ();
                @data = split /\t/, $dataline;
#               $col = @data;
#               next unless ($col == 5);
                next unless ($data[0] eq $chr);
                $bStart =  $data[2];
                $bStop = $data[3];

                if ($data[1] eq $strand) { #positive strand
#                       next unless ($stop > $bStart);
		        next unless ($lateStart > $bStart);
                	$dis = $start - $bStop;
                	next unless ($dis < 2000);			
			$foundTF = 1;
		} else { #negative strand
#                       next unless ($start < $bStop);
		        next unless ($lateStart < $bStop);
                        $dis = $bStart - $stop;
                        next unless ($dis < 2000);
                        $foundTF = 1;
		}
		next unless ($foundTF == 1);
                #found a TF for this gene
                $tf = $data[4];
                $gTF = join "", $g, $tf;
                if ($existTF{$gTF}) {
			#do nothing
                } else {
                   $existTF{$gTF} = 1;
                   if ($tfList eq "") {
                      $tfList = $tf;
                   } else {
                      $tfList = join ", ", $tfList, $tf;
                   }
               }
                       
        }
        close (DATA);
	if ($tfList eq "") {
	   $tfList = "N.A.";
	}
	my @allTF = split ", ", $tfList;
	@sortTF = ();
	@sortTF =  sort { lc($a) cmp lc($b) } @allTF;
	$tfList = join ", ", @sortTF;
	print OUT "$g\t$map\t$chr $strand $start $stop\t$tfList\n";	
    } else {
	#this gene is not on GFF3
	#print "$g is not on GFF3!\n";
        print OUT "$g\t$map\tN.A.\tN.A.\n";
    }
}

close (OUT);
print "Done printing the genetic map, chromosome coordinates and transcription factor table for $totalGenes genes.\n";
