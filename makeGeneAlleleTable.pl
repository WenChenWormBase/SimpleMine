#!/usr/bin/perl -w
use strict;

print "This script create Gene-Allele table based on an ace file dumped from the current WS release.\n"; # the dump only include sequenced allele variation objects

my %aGene; #var ID -> WB gene ID
my %aName; #var ID -> Public allele Name
my %aType; #var ID -> Mutation type (Substitution, Insertion, Deletion, Tandem_duplication)
my %aMolChange; #var ID -> Molecular_change (Missense, Nonsense, Splice_site, Frameshift)
my %aExon; #var ID -> Coding_exon, Intron, Promoter, UTR_3, UTR_5

my %geneAlleleList; #a list of alleles for each gene
my @geneWithAllele; #a list of genes that have sequenced alleles
my %existGene; #record if the gene is already documented in @geneWithAllele

my ($a, $g, $alleleName, $line);
my @tmp;

open (IN, "/home/wen/simpleMine/ace_files/SequencedAllele.ace") || die "can't open SequencedAllele.ace!";
#open (IN, "/home/wen/simpleMine/ace_files/testSequencedAllele.ace") || die "can't open SequencedAllele.ace!";


#----- print Gene-Allele table -----

my $i = 0; #for @geneAlleleList
while ($line = <IN>) {
    chomp ($line);
    if ($line =~ /^Variation/) {
	@tmp = ();
	@tmp = split '"', $line;
	$a = $tmp[1];
	#print "Processing allele ID $a.\n";
	$aGene{$a} = "N.A.";
	$aName{$a} = "N.A.";
	$aType{$a} = "N.A.";
	$aMolChange{$a} = "N.A.";
	$aExon{$a} = "N.A.";	
    } elsif ($line =~ /^Public_name/) {
	@tmp = ();
	@tmp = split '"', $line;
	$aName{$a} = $tmp[1]; #get the public allele name
    } elsif ($line =~ /^Gene/) {
	@tmp = ();
	@tmp = split '"', $line;
	$g = $tmp[1];
	$aGene{$a} = $g; #get the associated Gene ID
	if ($geneAlleleList{$g}) {
	    $geneAlleleList{$g} = join ",", $geneAlleleList{$g}, $a;
	} else {
	    $geneAlleleList{$g} = $a;
	}
	if ($existGene{$g}) {
	    #do nothing, this gene is already listed
	} else {
	    $geneWithAllele[$i] = $g;
	    $i++;
	    $existGene{$g} = 1;
	}	
    } elsif ($line =~ /^Substitution/) {
	$aType{$a} = "Substitution";
    } elsif ($line =~ /^Insertion/) {
	$aType{$a} = "Insertion";
    } elsif ($line =~ /^Deletion/) {
	$aType{$a} = "Deletion";
    } elsif ($line =~ /^Tandem_duplication/) {
	$aType{$a} = "Tandem_duplication";
    }
    
    if (($line =~ /^Predicted_CDS/)||($line =~ /^Gene/)||($line =~ /^Transcript/)||($line =~ /^Pseudogene/)) {
	if ($line =~ /Missense/) {
	    $aMolChange{$a} = "Missense";
	} elsif ($line =~ /Missense/) {
	    $aMolChange{$a} = "Missense";
	#} elsif ($line =~ /Nonsense/) {
	} elsif (($line =~ /Amber_UAG/)||($line =~ /Ochre_UAA/)||($line =~ /Opal_UGA/)||($line =~ /Ochre_UAA_or_Opal_UGA/)||($line =~ /Amber_UAG_or_Ochre_UAA/)||($line =~ /Amber_UAG_or_Opal_UGA/)) {    
	    $aMolChange{$a} = "Nonsense";
	} elsif ($line =~ /Splice_site/) {
	    $aMolChange{$a} = "Splice_site";
	} elsif ($line =~ /Frameshift/) {
	    $aMolChange{$a} = "Frameshift";
	} elsif ($line =~ /Silent/) {
	    $aMolChange{$a} = "Silent";
	    #print "$a is a silent allele.\n";
	} elsif ($line =~ /Coding_exon/) {
	    $aExon{$a} = "Coding_exon";
	#} elsif ($line =~ /Intron/) {
	#    $aExon{$a} = "Intron";
	#} elsif ($line =~ /Noncoding_exon/) {
	#    $aExon{$a} = "Noncoding_exon";
	#} elsif ($line =~ /Promoter/) {
	#    $aExon{$a} = "Promoter";
	#} elsif ($line =~ /UTR_3/) {
	#    $aExon{$a} = "UTR_3";
	#} elsif ($line =~ /UTR_5/) {
	#    $aExon{$a} = "UTR_5";
	}
    } 
}    

close (IN);
print "A total of $i genes have sequenced allele information.\n";


open (OUT, ">GeneAllele.csv") || die "cannot open $!\n";
#print OUT "WormBase Gene ID\tSequenced Allele\n";
print OUT "WormBase Gene ID\tCoding_exon Non_silent Allele";

#my @aList; # a list of alleles for each gene
#my $p = 0; #record if the allele is already printed
my @delAllele;
my @insAllele;
my @subAllele;
my @tanAllele;
my @sortDelA;
my @sortInsA;
my @sortSubA;
my @sortTanA;
my ($aDes, $delAi, $insAi, $subAi, $tanAi, $allAlleleDes);


foreach $g (@geneWithAllele) {

    #sort the list according to allele type
    @tmp = ();
#    @aList = ();
#    $i = 0; #this is for sorting @aList
    @tmp = split ",", $geneAlleleList{$g};

    @delAllele = ();
    @insAllele = ();
    @subAllele = ();
    @tanAllele = ();
    @sortDelA = ();
    @sortInsA = ();
    @sortSubA = ();
    @sortTanA = ();
    $delAi = 0;
    $insAi = 0;
    $subAi = 0;
    $tanAi = 0;
    
    foreach $a (@tmp) {
	next unless ($aName{$a});
	next unless ($aExon{$a} eq "Coding_exon");
	next unless ($aMolChange{$a} ne "Silent");
	
	$aDes = "$aName{$a}\|$aType{$a}\|$aMolChange{$a}";
	
	if ($aType{$a} eq "Deletion") {
	    $delAllele[$delAi] = $aDes;
	    $delAi++;
	} elsif ($aType{$a} eq "Insertion") {
	    $insAllele[$insAi] = $aDes;
	    $insAi++;
	} elsif  ($aType{$a} eq "Substitution") {
	    $subAllele[$subAi] = $aDes;
	    $subAi++;
	} elsif  ($aType{$a} eq "Tandem_duplication") {
	    $tanAllele[$tanAi] = $aDes;
	    $tanAi++;
	} 
    } 

    @sortDelA = sort @delAllele;
    @sortInsA = sort @insAllele;
    @sortSubA = sort @subAllele;
    @sortTanA = sort @tanAllele;

    $allAlleleDes = join ", ", @sortDelA, @sortInsA, @sortSubA, @sortTanA;
    print OUT "\n$g\t$allAlleleDes";
    
#    foreach $a (@tmp) {
#	next unless ($aName{$a});
#	next unless ($aExon{$a} eq "Coding_exon");
#	next unless ($aMolChange{$a} ne "Silent");
#	
#	next unless ($aType{$a} eq "Deletion");
#	$aList[$i] = $a;
#	$i++;
#    } 
#
#    foreach $a (@tmp) {
#	next unless ($aName{$a});
#	next unless ($aExon{$a} eq "Coding_exon");
#	next unless ($aMolChange{$a} ne "Silent");
#
#	next unless ($aType{$a} eq "Insertion");
#	$aList[$i] = $a;
#	$i++;
#    } 

#    foreach $a (@tmp) {
#	next unless ($aName{$a});
#	next unless ($aExon{$a} eq "Coding_exon");	
#	next unless ($aMolChange{$a} ne "Silent");
#
#	next unless ($aType{$a} eq "Substitution");
#	$aList[$i] = $a;
#	$i++;
#    }
    
#    foreach $a (@tmp) {
#	next unless ($aName{$a});
#	next unless ($aExon{$a} eq "Coding_exon");
#	next unless ($aMolChange{$a} ne "Silent");
#
#	next unless ($aType{$a} eq "Tandem_duplication");
#	$aList[$i] = $a;
#	$i++;
#   } 
#
    
    #print allele table

    
#    $p = 0;
#    foreach $a (@aList) {
#	if ($aMolChange{$a} ne "Silent") {
#	    $alleleName = $aName{$a};
#	    if ($p == 1) {
#		print OUT ", $alleleName\|$aType{$a}\|$aMolChange{$a}";
#	    } else {	
#		print OUT "\n$g\t$alleleName\|$aType{$a}\|$aMolChange{$a}";
#		$p = 1;
#	    }
#	}
#    }
}
close (OUT);
print "Done printing Gene-Allele table.\n";
