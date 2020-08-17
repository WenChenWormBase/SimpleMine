#!/usr/bin/perl -w

use strict;

print "Combine all tables and split them into species specific files.\n";

my ($line, $g, $f, $filename, $dataline, $naValue, $key, $totalFields, $fullspename, $pubname, $checkstring, $allothername, $allsynonym, $name);
my @tmp;
my @emptyValue;
my @otherNameList;
my @headerFields;
my %gNameExist; #check if a combination of gene id and name already exist
my %gName; #nameline for each gene entry;
my %gData; #dataline for each gene entry;
my @fileList = ("1_WBGeneName.csv", "2_GeneticMapPosition.csv", "3_RNAiAllelePhenoObserved.csv", "4_GeneAllele.csv", "5_ConfirmedInteraction.csv", "6_GeneTissueLifeStage.csv", "7_GeneDiseaseHumanOrtholog.csv", "8_GeneReference.csv", "9_1_GeneOntologyAssociation.csv", "9_2_GeneDescription.csv");
my $filepath1 = "/home/wen/simpleMine/sourceFile/";

my $i = 0; #for gene list
my $gapGene = 0; #for empty genes to fill gaps
my @genelist;
my %existGene;
my %foundValueGene; 
my $header = "";
my $nameheader;

foreach $f (@fileList) {
    $filename = join "", $filepath1, $f;
    print "Read $filename ...";
    
    open (IN, "$filename") || die "can't open $filename!";
    $line =<IN>;
    chomp ($line);
    @tmp = ();
    @tmp = split /\t/, $line;
    $key = shift @tmp;

    #creat N.A. entries
    @emptyValue = ();
    @emptyValue = @tmp;
    foreach (@emptyValue) {
	$_ = "N.A.";
    }
    $naValue = join "\t", @emptyValue;
    print "$naValue ...";
    
    if ($header eq "" ) {
	$header = join "\t", $key, @tmp;
	#$nameheader = $header;
	$nameheader = "WormBase Gene ID\tPublic Name\tSynonym";
    } else {
	$header = join "\t", $header, @tmp;
    }
    
    while ($line =<IN>) {
	chomp($line);
	@tmp = ();
	@tmp = split /\t/, $line;
	$g = shift @tmp;

	unless ($existGene{$g}) {
	    next unless ($f eq "1_WBGeneName.csv");
	    next unless !($line =~ /Trichuris muris/);
	    $genelist[$i] = $g;
	    $existGene{$g} = 1;
	    $i++;
	}
	
	if ($gData{$g}) {	    
	    $gData{$g} = join "\t", $gData{$g}, @tmp;
	    $foundValueGene{$g} = $f;
	} else {
	    
	    next unless ($f eq "1_WBGeneName.csv");
	    $gData{$g} = join "\t", @tmp;
	    $foundValueGene{$g} = $f;
	    
	    $pubname = shift @tmp;
	    $checkstring = join "", $g, $pubname;
	    $gNameExist{$checkstring} = 1;	    
	    $fullspename = shift @tmp;
	    my @syn = ();
	    my $si = 0;
	    @otherNameList = ();
	    $allothername = join ", ", @tmp;	    
	    @otherNameList = split ", ", $allothername; 
	    foreach $name (@otherNameList) {
		next unless ($name ne "");
		next unless ($name ne "N.A.");
		$checkstring = join "", $g, $name;
		next unless !($gNameExist{$checkstring});
		$syn[$si] = $name;
		$si++;
		$gNameExist{$checkstring} = 1;
	    }
	    $allsynonym = join ", ", @syn;
	    $gName{$g} = join "\t", $pubname, $allsynonym;
	    
	    #$gName{$g} = $gData{$g};
	}
    }    
    close (IN);
    print "done.\n";
    
    #fill the gap with N.A. fields	    
    next unless ($f ne "1_WBGeneName.csv");
    $gapGene = 0;
    foreach $g (@genelist) {
	if ($foundValueGene{$g}){
	    if ($foundValueGene{$g} ne $f) {
		$gData{$g} = join "\t", $gData{$g}, $naValue;
		$foundValueGene{$g} = $f;
		$gapGene++;
	    }
	} else {
	    print "$f - no data for $g\n.";
	    #$gData{$g} = $naValue;
	    #$foundValueGene{$g} = $f;
	    #$gapGene++;
	}
    }
    print "fill empty gaps for $gapGene genes.\n";
    
}
print "Processed a total of $i genes.\n";
#print "$header\n";


#---- print mixed species table ---
my $filepath2 = "/home/wen/simpleMine/multiSpeSimpleMine/";
#my @speDir = ("b_malayi", "c_brenneri", "c_briggsae", "c_elegans", "c_japonica", "C_remanei", "o_volvulus", "p_pacificus", "s_ratti", "mix");
my @speDir = ("b_malayi", "c_brenneri", "c_briggsae", "c_elegans", "c_japonica", "c_remanei", "o_volvulus", "p_pacificus", "s_ratti");
my %fullSpeName = ("b_malayi" => "Brugia malayi", 
		   "c_brenneri" => "Caenorhabditis brenneri", 
		   "c_briggsae" => "Caenorhabditis briggsae", 
		   "c_elegans" => "Caenorhabditis elegans", 
		   "c_japonica" => "Caenorhabditis japonica", 
		   "c_remanei" => "Caenorhabditis remanei", 
		   "o_volvulus" => "Onchocerca volvulus", 
		   "p_pacificus" => "Pristionchus pacificus", 
		   "s_ratti" => "Strongyloides ratti");


my $mixnamefile = join "", $filepath2, "mix", "\/", "GeneName.csv";
my $mixdatafile = join "", $filepath2, "mix", "\/", "SimpleMineSourceData.csv"; 
open (MIXNAM, ">$mixnamefile") || die "cannot open $mixnamefile!\n";
print MIXNAM "$nameheader\n";
open (MIXDAT, ">$mixdatafile") || die "cannot open $mixnamefile!\n";
print MIXDAT "$header\n";

foreach $f (@speDir) {

    my $datafile = join "", $filepath2, $f, "\/", "SimpleMineSourceData.csv";
    my $namefile = join "", $filepath2, $f, "\/", "GeneName.csv";
    
    open (NAM, ">$namefile") || die "cannot open $namefile!\n";
    print NAM "$nameheader\n";
    
    open (DAT, ">$datafile") || die "cannot open $namefile!\n";
    print DAT "$header\n";

    $i = 0;
    next unless $fullSpeName{$f};
    $fullspename = $fullSpeName{$f};
    
    print "work on $fullspename ...";
    
    foreach $g (@genelist) {
	if ($gData{$g}) {	    
	    
	    $dataline = $gData{$g};
	    @tmp = ();
	    @tmp = split /\t/, $dataline;
	    $totalFields = @tmp;
	    next unless ($totalFields > 2);
	    my $spe = $tmp[1];
	    next unless ($spe eq $fullspename);
	    print NAM "$g\t$gName{$g}\n";
	    print DAT "$g\t$gData{$g}\n";	    
	    print MIXNAM "$g\t$gName{$g}\n";
	    print MIXDAT "$g\t$gData{$g}\n";
	    $i++;

	}
	
    }
    
    print "$i genes parsed.\n";
    close (DAT);
    close (NAM);
}
close (MIXNAM);
close (MIXDAT);

my $headerfile = join "", $filepath2, "headers";
open (OUT, ">$headerfile") || die "cannot open $headerfile!\n";
print OUT "$header\n";
close (OUT);
