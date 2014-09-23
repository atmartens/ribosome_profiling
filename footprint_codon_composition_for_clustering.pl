#!/usr/bin/perl -w
use strict;

# Compute the codon occurrences of ribosome footprints overall, dependent both on
# position in the footprint and categorized by the length of the footprint
# 
# Use as input a "simplified" version of a SAM file
# 1. sequences 2. mapped position 3. chromosome id 4. number of times
# 
# ignore lengths of footprints, align all footprints to the 3' end (Right Hand Side)
# remove the leftmost codon because of nuclease signal
# only works with single chromosome files
# 
# Output file is formatted to work with M. de Hoon's Cluster program,
# but if interested in the length-independent, 3' aligned codon frequencies
# can just look at a particular column for a given position

my $fnaFile = $ARGV[0];
my $pttFile = $ARGV[1];
my $fpsFile = $ARGV[2];
my $outFile = $ARGV[3];

# Store the whole genome file which is assumed to be fna format
open (INGENOME, "<$fnaFile") or die $!;

# discard the first line, assumed to be descriptive
<INGENOME>;

my $genome = "";
while (<INGENOME>) {
    chomp;
    $genome .= $_;
}
close INGENOME;

# file containing gene coordinates + names: use a PTT file
open (COORDINATES, "<$pttFile") or die $!;
my @coords = <COORDINATES>;
close COORDINATES;

# drop the first 3 lines of the PTT file
for (0 .. 2) {
    shift @coords;
}

# store the frame information; initialize to zero
my @frames;
for (0 .. length($genome)) {
    $frames[$_] = 0;
}

# annotate across the length of the chromosome,
# for each position, whether it is in the 1st, 2nd, or 3rd position
# in the codon (aka the frame)
# the entire genome of the organism
foreach (@coords) {
    # location, strand, length, PID, gene, synonym, code, COG, product
    my @line = split /\t/;
    
    my $name = $line[4];
    my $prod = $line[8];
    
    my $range = $line[0];
    my $strand = $line[1];
    
    # parse the range to get the beginning and end
    my ($left, $right) = split /\.\./, $range;
    
    if ($strand eq "+") {
	for (my $i = $left; $i < $right; $i += 3) {
	    $frames[$i] = 1;
	    $frames[$i+1] = 2;
	    $frames[$i+2] = 3;
	}
    }
    
    else {
	for (my $i = $left; $i < $right; $i += 3) {
	    $frames[$i] = -3;
	    $frames[$i+1] = -2;
	    $frames[$i+2] = -1;
	}
    }
}

# All 64 codons
my @codonlist = @{CodonList()};

# the number of occurrences
# store a list for each position
my @counts;

# footprints file in simplified format
open (IN, "<$fpsFile") or die $!;

LINE: while (<IN>) {
    # sequence, mapped position, chromosome ID, counts
    # need to determine if frame is correct
    chomp;
    my @line = split /\t/;
    
    my $sequence = $line[0];
    my $startpos = $line[1];
    my $num = $line[3];
    my $l = length($sequence);
    
    my $seq_inframe;
    
    # check every single base to make sure no frame equals zero
    foreach my $base ($startpos .. $startpos + $l) {
	if ($frames[$base] == 0) {
	    # skip the whole line in the while loop
	    next LINE;
	}
    }
    
    # this frame-detection strategy only removes from the beginning and ignores the end
    # this works because we only care about triplets later on
    
    # positive strand
    if ($frames[$startpos] > 0) {
	my $remove = (4 - $frames[$startpos]) % 3;
	$seq_inframe = substr $sequence, $remove;
    }
    
    # negative strand
    else { 
	my $remove = (4 + $frames[$startpos + $l - 1]) % 3;
	$seq_inframe = substr $sequence, $remove;
    }
    
    # convert sequence to list of triplets
    my @codons = $seq_inframe =~ /([ATGC]{3})/g;
    
    # remove the leftmost codon
    # otherwise we are mixing in the bias at the 5' end for the diff. lengths!
    shift @codons;
    
    # reverse because that way they are all indexed at zero
    # and it makes it easier to normalize by position
    @codons = reverse @codons;
    
    # Count the codons
    # scale by the number of times that this footprint was found ~ aka in $num
    for my $pos (0 .. $#codons) {
	unless (exists($counts[$pos])) {
	    my %newhash;
	    $counts[$pos] = \%newhash;
	}
	$counts[$pos]{$codons[$pos]} += $num;
    }
}

close IN;

# Need to go through the tabulations and compute,
# at each position, the total number of counts across all codons.
# Loop through the existing data structures and store these values in
# new data structures

# store the total number codons across all positions
my @sums;

foreach my $pos (0 .. $#counts) {
    foreach my $c (@codonlist) {
	if (defined($sums[$pos])) {
	    $sums[$pos] += $counts[$pos]{$c};
	}
	# initialize if doesn't yet exist
	else {
	    $sums[$pos] = $counts[$pos]{$c};
	}
    }
}

open (OUT, ">$outFile") or die $!;
# some clustering/treeviewing programs demand YORF\tNAME
print OUT "NAME\tNAME";
for (my $j = scalar(@counts) - 1; $j > 0; $j--) {
    print OUT "\t$j";
}
print OUT "\n";

foreach my $c (@codonlist) {
    print OUT "$c\t$c";
    # minus one more because lack of data in the final column
    for (my $pos = scalar(@counts) - 2; $pos >= 0; $pos--) {
	print OUT "\t";
	# calculate the percentage by dividing the counts at a position
	# by the total at a position, then multiply by 100
	print OUT ($counts[$pos]{$c} / $sums[$pos] * 100);
    }
    print OUT "\n";
}
close OUT;


# Sub #

# make all codons
# returns a reference to the newly-created list
# because the sequencing data are in DNA, make these codons in DNA

sub CodonList {
    my @letters = ('T','C','A','G');
    my @codon_list;

    for (my $i = 0; $i < 4; $i++) {
	for (my $j = 0; $j < 4; $j++) {
	    for (my $k = 0; $k < 4; $k++) {
		my $c = $letters[$i] . $letters[$j] . $letters[$k];
		  push @codon_list, $c;
	    }
	}
    }
    
    return \@codon_list;
}
