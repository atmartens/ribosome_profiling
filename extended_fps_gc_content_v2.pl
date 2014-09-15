#!/usr/bin/perl -w
use strict;

# Compute GC content counts from an "extended" set of footprints
# assume all seqs are extended to length 45
# take as input a file with:
# 1. extended sequence
# 2. chromosome ID
# 3. number of occurrences
# 4. original length

my @nucleotides = ("A", "T", "G", "C");

my $inFile = $ARGV[0];

# Read the sequences and tally the nucleotides per position
# by length category (original footprint length before extension)

open (IN, "<$inFile") or die $!;

my %lengthCategory;

while (<IN>) {
    my @line = split;
    my $s = $line[0];
    my $n = $line[2];
    my $l = $line[3];
    
    # Create hash of ATGC counts if it doesn't exist already
    # for this length category
    unless (defined($lengthCategory{$l})) {
	my @positions;
	$lengthCategory{$l} = \@positions;
	for my $p (0 .. 44) {
	    my %counts;
	    foreach my $nuc (@nucleotides) {
		$counts{$nuc} = 0;
	    }
	    $positions[$p] = \%counts;
	}
    }
    
    for my $p (0 .. 44) {
	my $base = substr $s, $p, 1;
	$lengthCategory{$l}[$p]{$base} += $n;
    }
}

close IN;

# Go from smallest to greatest length
my @high_low_num = sort {$a <=> $b} keys %lengthCategory;

open (OUT, ">$ARGV[1].txt") or die $!;

# Print out lines for 0 to 20 which are empty
# so that the Y-axis in the heatmap starts at 21
for (0 .. 20) {
    for (0 .. 43) {
        print OUT "-\t";
    }
    print OUT "-\n";
}

for my $length (@high_low_num) {
    # total number of bases - assume is same at all positions
    my $total = 0;
    
    foreach my $letter (keys $lengthCategory{$length}[0]) {
	$total += $lengthCategory{$length}[0]{$letter};
    }
    
    for my $p (0 .. 44) {
	my $GC = ($lengthCategory{$length}[$p]{"G"} + $lengthCategory{$length}[$p]{"C"}) / $total;
	print OUT "$GC\t";
    }
    
    print OUT "\n";
}

open GNUPLOT, "| gnuplot";
print GNUPLOT <<gnuplot_commands;
set terminal png
set output "$ARGV[1].png"
# set title "GC content in ribosomal footprints of different lengths"
set xlabel "Position in footprint (5' to 3')"
set ylabel "Footprint length"
set xrange [-0.5:44.5]
set yrange [20.5:45.5]
set xtic nomirror
set ytic nomirror
plot "$ARGV[1].txt" matrix with image
unset output
gnuplot_commands
close GNUPLOT;
