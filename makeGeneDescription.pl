#!/usr/bin/perl -w

use strict;

print "This script create Gene-Desciption table based on an ace file dumped from the current WS release.\n";

my $line;
my $g;
my $ecs;
my $autoDes;
my $conDes;
my @tmp;
my $i = 0;
my $totalGenes;

#------- Get Expression Cluster Summary ---------------------------------
open (EC, "/home/wen/AutoDescription/ecSummary/allECreg.csv") || die "can't open allECreg.csv!";
my %ecSum;
while ($line = <EC>) {
    next unless ($line =~ /^WBGene/);
    chomp ($line);
    @tmp = split /\t/, $line;
    $g = $tmp[0];
    $ecs = $tmp[2];
    if ($ecSum{$g}) {
	$ecSum{$g} = join " ", $ecSum{$g}, $ecs;
    } else {
	$ecSum{$g} = $ecs;
    }
}
close (EC);

#----------------------- Get Automated and Concise Descriptions from WS ---------------
open (IN, "/home/wen/simpleMine/ace_files/WBGeneDescription.ace") || die "can't open WBGeneDescription.ace!";
open (OUT, ">GeneDescription.csv") || die "cannot open $!\n";
print OUT "WormBase Gene ID\tConcise Description\tAutomated Description\tExpression Cluster Summary\n";

while ($line = <IN>) {
    chomp ($line);
    next unless ($line ne "");
    @tmp = ();
    if ($line=~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	$autoDes = "N.A.";
	$conDes = "N.A.";
	$i++;
    } elsif ($line =~ /^Concise_description/) {
	@tmp = split '"', $line;
	$conDes = $tmp[1];	
    } elsif (($line =~ /^Automated_description/)&&($line =~ /Date_last_updated/)) {

	@tmp = split '"', $line;
	$autoDes = $tmp[1];

	#---- get EC Summary ----
	if ($ecSum{$g}) {
	    $ecs = $ecSum{$g};
	} else {
	    $ecs = "N.A.";
	}

	#---- done getting EC Summary ---
	print OUT "$g\t$conDes\t$autoDes\t$ecs\n";
    }
}    
$totalGenes = $i;
close (IN);
close (OUT);
print "Done printing Gene-Description table for $totalGenes genes.\n";
