#!/usr/bin/perl -w

use strict;
use Ace;

print "This script create Gene-RNAi-Allele Phenotype table from current WS release.\n";

my $acedbpath='/home/citace/WS/acedb/';
my $tace='/usr/local/bin/tace';

print "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) || die print "Connection failure: ", Ace->error;

#------------------ find all genes that contains RNAi or Allele ----
my ($spe, $ref, $t, $i, $r, $pname, $all_gene, $all_phenotype, $no_phenotype, $g);
my @tmp;
my $totalResult;
my %gRNAi;
my %gAllele;

my $query="QUERY FIND Gene RNAi_result = *";
my @rnaiGene = $db->find($query);
foreach $g (@rnaiGene) {
    $gRNAi{$g} = 1;
}
print scalar @rnaiGene, " genes found with RNAi results.\n";

$query="QUERY FIND Gene Allele = *";
my @alleleGene = $db->find($query);
foreach $g (@alleleGene) {
    $gAllele{$g} = 1;
}
print scalar @alleleGene, " genes found with Allele.\n";

#------------------Build Phenotype Description table ------------

open (OUT, ">RNAiAllelePheno.csv") || die "cannot open $!\n";
print OUT "Gene\tRNAi Phenotype Observed\tRNAi Phenotype Not Observed\tAllele Phenotype Observed\tAllele Phenotype Not Observed\n";

$query="QUERY FIND Gene WBGene*";
my @gene = $db->find($query);
print scalar @gene, " genes total found in ACeDB.\n";

my ($RNAiPhe, $RNAiNotPhe, $AllelePhe, $AlleleNotPhe);

foreach $g (@gene) {
    if ($gRNAi{$g}) {
	#get RNAi Phenotype
	$query="QUERY FIND Gene $g; follow RNAi_result; follow Phenotype; follow Primary_name"; 
	@tmp = ();
	@tmp = $db->find($query);
	$totalResult = @tmp;
	if ($totalResult == 0) {
	    $RNAiPhe = "N.A.";
	} else {
	    $RNAiPhe = join ",", @tmp;
	}

	#get RNAi Not_observed phenotype
	$query="QUERY FIND Gene $g; follow RNAi_result; follow Phenotype_not_observed; follow Primary_name"; 
	@tmp = ();
	@tmp = $db->find($query);
	$totalResult = @tmp;
	if ($totalResult == 0) {
	    $RNAiNotPhe = "N.A.";
	} else {
	    $RNAiNotPhe = join ",", @tmp;
	}
    } else {
	$RNAiPhe = "N.A.";
	$RNAiNotPhe = "N.A.";
    }

    if ($gAllele{$g}) {
	#get Allele Phenotype
	$query="QUERY FIND Gene $g; follow Allele; follow Phenotype; follow Primary_name"; 
	@tmp = ();
	@tmp = $db->find($query);
	$totalResult = @tmp;
	if ($totalResult == 0) {
	    $AllelePhe = "N.A.";
	} else {
	    $AllelePhe = join ",", @tmp;
	}

	#get Allele Not_observed Phenotype
	$query="QUERY FIND Gene $g; follow Allele; follow Phenotype_not_observed; follow Primary_name"; 
	@tmp = ();
	@tmp = $db->find($query);
	$totalResult = @tmp;
	if ($totalResult == 0) {
	    $AlleleNotPhe = "N.A.";
	} else {
	    $AlleleNotPhe = join ",", @tmp;
	}
    } else {
	$AllelePhe = "N.A.";
	$AlleleNotPhe = "N.A.";
    }
    print OUT "$g\t$RNAiPhe\t$RNAiNotPhe\t$AllelePhe\t$AlleleNotPhe\n";
    
}

close (OUT);
$db->close();
