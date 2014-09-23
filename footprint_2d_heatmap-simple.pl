#!/usr/bin/perl -w
use strict;

# Take footprinting data and produce a heatmap representing the 2 ends of the footprints
# along the gene by heatmap of counts
# be sure to remove the output from a previous run, since we are APPENDING
# This is part 1: generates files containing the beginning, ending positions for each read
# for each gene.

# NOTE: This ignores which strand a gene is on, and only tallies
# the absolute chromosome positions all in the same direction from the zero point
# for now, just handle a single pair of FNA/PTT files
# NOTE: issues with overlapping ORFs?

my $inGenome = $ARGV[0];
my $inPTT = $ARGV[1];
my $inReads = $ARGV[2];
my $outDir = $ARGV[3];

# the entire genome of the organism
# Store the whole genome file which is assumed to be fna format
open (INGENOME, "<$inGenome") or die $!;

# discard the first line, assumed to be descriptive
<INGENOME>;

my $genome = "";
while (<INGENOME>) {
    chomp;
    $genome .= $_;
}
close INGENOME;


# file containing gene coordinates + names: use a PTT file
open (COORDINATES, "<$inPTT") or die $!;
my @coords = <COORDINATES>;
close COORDINATES;

# drop the first 3 lines of the PTT file
for (0 .. 2) {
    shift @coords;
}

my @gene_names;
my %endpoints;

foreach (@coords) {
    chomp;
    
    my @line = split /\t/;
    # location, strand, length, PID, gene, synonym, code, COG, product
    
    my $name = $line[4];
    my $code = $line[5];
    my $prod = $line[8];
    
    my $range = $line[0];
    my $strand = $line[1];
    
    # parse the range to get the beginning and end
    my ($left, $right) = split /\.\./, $range;
    
    my $name_code = "$name-$code";
    
    # also use the code, since this is unique; some of the names
    # are repeated for diff. genes, such as transposon insertions
    push @gene_names, $name_code;
    my @temparr = ($left, $right);
    $endpoints{$name_code} = \@temparr;
}

# Go through the *simplified* SAM file, and for each footprint,
# extract its beginning, end, and store this in a file
# for the corresponding gene

open (SAM, "<$inReads") or die $!;

my $outdir = $outDir;

while (<SAM>) {
    my @line = split;
    # $line[0] : the sequence of the read
    # $line[1] : the chromosomal position of the alignment
    # $line[2] : name of chromosome
    # $line[3] : num times occurs
    
    my $read = $line[0];
    my $startpos = $line[1];
    my $endpos = $startpos + length($read);
    my $num = $line[3];
    
    
    # what to do if the footprint goes over the end of the gene? for now, ignore it
    foreach my $n (@gene_names) {
	if (($startpos >= ${$endpoints{$n}}[0]) && ($endpos <= ${$endpoints{$n}}[1])) {
	    open (OUT, ">>$outdir/$n.txt") or die $!;
	    for (0 .. $num) {
		print OUT "$startpos\t$endpos\n";
	    }
	    close OUT;
	    
	    # this will short-circuit the foreach loop once we find the right gene
	    # with the starting and ending points in the correct range
	    last;
	}
    }
}

close SAM;
