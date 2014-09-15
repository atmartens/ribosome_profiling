#!/usr/bin/perl -w
use strict;

# Take a directory of "endpoint" files which are series of number denoting the starting and 
# ending genomic coordinates of footprints for a single gene.
# Collapse them into files which have the number of each pair of coordinates.

my $inDir = $ARGV[0];
my $outDir = $ARGV[1];

my @inFiles = `ls -1 $inDir`;

foreach my $file (@inFiles) {
    chomp $file;
    my %coords;
    
    open(IN, "<$inDir/$file") or die $!;
    while (<IN>) {
        chomp;
        if (defined($coords{$_})) {
            $coords{$_}++;
        }
        else {
            $coords{$_} = 1;
        }
    }
    close IN;
    
    open (OUT, ">$outDir/$file") or die $!;
    foreach (keys %coords) {
        print OUT "$_\t$coords{$_}\n";
    }
    close OUT;
}
