#!/usr/bin/perl -w

use strict;
use Ace;

my @AO;
my @LS;
#my %geneName;
my ($g, $totalTissue, $totalGene, $totalLifeStage);

print "This script create Gene-Tissue-LifeStage table based on Expr_pattern data from current WS release.\n";

my $acedbpath='/home/citace/WS/acedb/';
my $tace='/usr/local/bin/tace';

print "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) || die print "Connection failure: ", Ace->error;

#------------Build expression cluster and expression pattern gene hashes-------
my $query="QUERY FIND Expression_cluster Anatomy_term = * AND NEXT = Enriched; follow Gene";
my @ecGene = $db->find($query);
my %gEC;
foreach $g (@ecGene) {
    $gEC{$g} = 1;
}
print scalar @ecGene, " genes contain anatomy enriched info from expression cluster.\n";

$query="QUERY FIND Expr_pattern Expr*; !Type = Microarray; !Type = EPIC; follow Gene;";
my @epGene = $db->find($query);
my %gEP;
foreach $g (@epGene) {
    $gEP{$g} = 1;
}
print scalar @epGene, " genes contain anatomy info from expression pattern.\n";

$query="QUERY FIND Gene WBGene*";
my @gene = $db->find($query);
print scalar @gene, " genes total found in ACeDB.\n";
#--------------Done ---------------------------------------------------


#-----------------Build Gene - Tissue table---------------------------
open (OUT1, ">GeneTissueLifeStage.csv") || die "cannot open $!\n";
print OUT1 "WormBase Gene ID\tExpr_pattern Tissue\tGenomic Study Tissue\tExpr_pattern LifeStage\tGenomic Study LifeStage\n";

my ($epAOlist, $epLSlist, $ecAOlist, $ecLSlist);
foreach $g (@gene) {

    #---- get anatomy_term from Expr_pattern ------------------
    if ($gEP{$g}) {
	@AO = ();
	$query="find Gene $g; follow Expr_pattern; !Type = Microarray; !Type = EPIC; follow Anatomy_term; follow Term";
	@AO = $db->find($query);
	$totalTissue = @AO;
 
	if  ($totalTissue == 0) {
	    $epAOlist = "N.A.";
	} else {
	    $epAOlist = join ", ", @AO; 
	}
    } else {
	$epAOlist = "N.A.";
    }
    

    #----- get anatomy_term from Expression_cluster ---------
    if ($gEC{$g}) {
	@AO = ();
	$query="find Gene $g; follow Expression_cluster; Anatomy_term = * AND NEXT = Enriched; follow Anatomy_term; follow Term";
	@AO = $db->find($query);
	$totalTissue = @AO;
 
	if  ($totalTissue == 0) {
	    $ecAOlist = "N.A.";
	} else {
	    $ecAOlist = join ", ", @AO; 
	}
    } else {
	$ecAOlist = "N.A.";
    }

    #----------- get life_stage from Expr_pattern ----------
    if ($gEP{$g}) {
	@LS = ();
	$query="find Gene $g; follow Expr_pattern; !Type = Microarray; !Type = EPIC; follow Life_stage; follow Public_name";
	@LS = $db->find($query);
	$totalLifeStage = @LS;
 
	if  ($totalLifeStage == 0) {
	    $epLSlist = "N.A.";
	} else {
	    $epLSlist = join ", ", @LS;
	    #print OUT1 "$LSlist";
	}
    } else {
	$epLSlist = "N.A.";
    }


    #----------- get life_stage from Expression_cluster ----------
     if ($gEC{$g}) {
	 @LS = ();
	 $query="find Gene $g; follow Expression_cluster; follow Life_stage; follow Public_name";
	 @LS = $db->find($query);
	 $totalLifeStage = @LS;
 
	 if  ($totalLifeStage == 0) {
	     $ecLSlist = "N.A.";
	 } else {
	     $ecLSlist = join ", ", @LS;
	 }
     } else {
	 $ecLSlist = "N.A.";
     }

    print OUT1 "$g\t$epAOlist\t$ecAOlist\t$epLSlist\t$ecLSlist\n";
}
close (OUT1);
print "Done printing Gene-Tissue-LifeStage table.\n";

$db->close();
