#!/usr/bin/perl -w

use strict;

print "This script create Gene-Map table based on an ace file dumped from the current WS release.\n";

my $line;
my $g;
my $chr;
my $p;
my @tmp;
my @pos;

open (IN, "/home/wen/simpleMine/ace_files/WBGeneMap.ace") || die "can't open WBGeneMap.ace!";


open (OUT, ">GeneticMapPosition.csv") || die "cannot open $!\n";
print OUT "WormBase Gene ID\tGenetic Map Position\n";

while ($line = <IN>) {
    chomp ($line);
    next unless ($line ne "");

    @tmp = ();
    @pos = ();
    if ($line=~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
    } elsif ($line =~ /^Map/) {
	@tmp = split '"', $line;
	$chr = $tmp[1];

	if ($line =~ /Position/) {
	    @pos = split /\s+/, $line;
	    $p = @pos[3];
	} else {
	    $p = "";
	}

	if ($p eq "") {
	    print OUT "$g\t$chr\n";
	} else {
	    print OUT "$g\t$chr $p\n";
	}
    }   
}    

close (IN);
close (OUT);
print "Done printing Gene-Map table.\n";
