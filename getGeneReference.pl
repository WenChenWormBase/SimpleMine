#!/usr/bin/perl -w

use strict;
use Ace;

print "This script create Gene-Reference table from the current WS release.\n";

my $acedbpath='/home/citace/WS/acedb/';
my $tace='/usr/local/bin/tace';

print "connecting to database... \n";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) || die print "Connection failure: ", Ace->error;


#----------- Build PMID - WBID hash ----------------------
my $query="query find Paper Database = MEDLINE";
my @pmidPaperList=$db->find($query);

my %PaperID;
my ($totalPaper, $pmid, $paper, $database);
$totalPaper = @pmidPaperList;

foreach $paper (@pmidPaperList) {
        $database = $paper->get('Database', 2);
	if ($database eq "PMID") {
	    $pmid=$paper->get('Database', 3);
	    $PaperID{$paper} = $pmid;
	} 
}
print "$totalPaper WormBase papers found with medline accession number.\n"; 


#------------Build disease term and id hash-------

#my %geneRef;
#my %geneRefCount;
my %paperPrimary;
my ($gene, $ref, $totalRef, $geneRef, $genePubMedRef);

#get the list of all primary research articles
$query="QUERY FIND Paper Type = Journal_article";
my @paperList = $db->find($query);
foreach $paper (@paperList) {
    $paperPrimary{$paper} = "PRA"; #primary research article
}


#get the list of genes associated with primary research articles
open (OUT, ">GeneReference.csv") || die "cannot open $!\n";
print OUT "WormBase Gene ID\tReference Count\tWB Paper ID\tPubMed ID\n";

$query="QUERY FIND Paper Type = Journal_article; follow Gene";
my @geneList = $db->find($query);
foreach $gene (@geneList) {
    $totalRef = 0;
    $geneRef = "";
    my @geneRefList = $gene -> Reference; 
    foreach $ref (@geneRefList) {
	next unless ($paperPrimary{$ref});

	if ($PaperID{$ref}) {
	  $pmid = $PaperID{$ref};
	} else {
	  $pmid = "N.A.";
	}
	
        if ($geneRef ne "") {
	   $geneRef = join ",", $geneRef, $ref;
	   $genePubMedRef = join ",", $genePubMedRef, $pmid;	    
	} else {
	    $geneRef = $ref;	   
	    $genePubMedRef = $pmid;
	}
        $totalRef++;
	
    }    
    print OUT "$gene\t$totalRef\t$geneRef\t$genePubMedRef\n";
}


close (OUT);
print "Done printing Gene-Reference table.\n";
$db->close();

