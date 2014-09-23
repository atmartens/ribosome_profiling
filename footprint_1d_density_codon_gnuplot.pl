#!/usr/bin/perl -w
use strict;

# Take tabulations of left to right footprints for a single gene
# Compute total occupancies of the codons, aka nucleotides divided by 3
# This is the smooshing of the 3' and 5' ends heatmaps plot to a 1D plot
# Make plots.

# the initial input files are series of start and endpoints of footprints
my @inputs = `ls -1 $ARGV[0]`;

foreach my $file (@inputs) {
    chomp $file;
    my $openfile = $ARGV[0] . $file;
    open (FILE, "<$openfile") or die $!;
    my @fps = <FILE>;
    close FILE;
    
    my %density;
    
    # tabulate
    foreach (@fps) {
	my ($lhs,$rhs) = split;
	for my $pos ($lhs .. $rhs) {
	    my $codon = int($pos / 3);
	    $density{$codon}++;
	}
    }
    
    # output
    my $extension = index($file, ".txt");
    my $name = substr $file, 0, $extension;
    my $outfile = "$ARGV[1]/$name-density.txt";
      
    open (OUT,">$outfile") or die $!;
    
    my @k = sort {$a <=> $b} keys %density;
    my $lowest = $k[0] - 1;
    
    foreach (@k) {
	print OUT $_ - $lowest . "\t$density{$_}\n";
    }
    
    close OUT;
    
    # run gnuplot on the output
    my $file_in = $outfile;
    my $file_out = "$ARGV[1]/$name-density.png";
    
    open GNUPLOT, "| gnuplot";
    print GNUPLOT <<gnuplot_commands;
    set terminal png
    set output "$file_out"
    set key off
    set title "$name Footprint Density"
    set ylabel "Number of Footprints"
    set xlabel "Codon Number"
    set autoscale
    set border 11
    set xtic nomirror
    plot "$file_in" using 1:2 with lines
    unset output
gnuplot_commands
    close GNUPLOT;
}
