#!/usr/bin/perl -w

use strict;

print "This script create Gene-GO_term table based on gene_association.WSXXX.wb and gene_ontology.WSXXX.obo files from the current WS release.\n";

my ($line, $g, $go, $goName, $checkstring, $goDes);
my @tmp;
my %goTerm; # names for each GO_term
my %geneGO; #a list of GO_terms for each gene.
my %existGeneGO; #this record if the GO_term is already documented for this gene
my @geneList; # all the genes that have go term
my $i;

#---- get GO names for all GO_terms ------------
open (IN1, "/home/wen/simpleMine/ace_files/gene_ontology.obo") || die "can't open gene_ontology.obo!";
while ($line =<IN1>) {
    chomp($line);
    if ($line =~ /^id/) {
	@tmp = split ": ", $line;
	$go = $tmp[1];
    } elsif ($line =~ /^name: /) {
	@tmp = split ": ", $line;
	$goName = $tmp[1];
	$goTerm{$go} = "\"$goName\"";
    }
}
close (IN1);


#------ go through gene association file ------------
open (IN2, "/home/wen/simpleMine/ace_files/gene_association.wb")  || die "cannot open gene_association.wb!\n";
$i = 0;
while ($line=<IN2>) {
    chomp($line);
    next unless ($line =~ /^WB/);
    @tmp = ();    
    @tmp = split /\t/, $line;
    $g = $tmp[1];
    $go = $tmp[4];
    $checkstring = join "", $g, $go;
    next unless !($existGeneGO{$checkstring}); #skip GO_term that was already recorded 
    next unless ($goTerm{$go});
    $existGeneGO{$checkstring} = 1;
    #$goDes = join ":", $go, $goTerm{$go};
    $goDes = $goTerm{$go};    
    if ($geneGO{$g}) {
	$geneGO{$g} = join "-----", $geneGO{$g}, $goDes; 
    } else {
	$geneGO{$g} = $goDes;
	$geneList[$i] = $g;
	$i++;
    }
}
close (IN2);
print "$i genes found with GO_term.\n";

#---- print Gene-GO_term table ----
open (OUT, ">GeneOntologyAssociation.csv") || die "cannot open GeneOntologyAssociation.csv!\n";
print OUT "WormBase Gene ID\tGene Ontology Association\n";

my $allGoNames;
my @goNameList;
my @sortedGo;
foreach $g (@geneList) {
    if ($geneGO{$g}) {
	@goNameList = ();
	@sortedGo = ();
	@goNameList = split "-----", $geneGO{$g};
	@sortedGo = sort { lc($a) cmp lc($b) } @goNameList;
	$allGoNames = join ", ", @sortedGo;
	
	print OUT "$g\t$allGoNames\n";
    } else {
	print "ERROR: $g has no ontology association!\n"
    }
}
close (OUT);
print "Done printing Gene-GO_term table.\n";

