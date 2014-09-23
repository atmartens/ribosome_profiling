#!/usr/bin/perl -w
use strict;

# Given a directory containing codon position-specific tallies of mRNA sequencing read counts
# for many genes, compute the average over the whole gene.
# Write all the averages to a new file.

# ARGUMENTS:
# 0. PTT file
# 1. FFN file
# 2. mRNA seq directory
# 3. output file for rowstacked histogram
# 4. output file for each codon alone
# 5. output file with ID and avg value

# Establish data structures:
# 1. Go through the PTT file and make a hash which inputs a unique identifier (e.g. b1276)
# and returns a genomic range
# 2. Go through the genome ffn file and make a hash which inputs a genomic range
# and outputs an mRNA sequence
# Next:
# 3. compute the average mRNA seq values and do the weighted codon counts

# UPDATE: remove first, last 10 codons from genes, because I do that for footprints too

##### 1.
# Read in the whole PTT file into a list
my $pttFile = $ARGV[0];
open (PTT, "<$pttFile") or die $!;
my @ptt = <PTT>;
close PTT;

# Remove the first 3 lines
for (0 .. 2) {
shift @ptt;
}

my %IDtoRange;

foreach (@ptt) {
    my @line = split;
    my $synonym = $line[5];
    my $range = $line[0];
    my ($left, $right) = split /\.\./, $range;
    my $strand = $line[1];
    
    my $loc;
    # If on minus strand, need to flip-flop and add a c
    if ($strand eq "-") {
	my $temp = $right;
	$right = $left;
	$left = $temp;
	$loc = "c" . $left . "-" . $right;
    }
    else {
	$loc = "$left-$right"
    }
    
    $IDtoRange{$synonym} = $loc;
}

##### 2. 
my $geneFile = $ARGV[1];

open (GENES, "<$geneFile") or die $!;

my %RangeToSeq;

my $r = "";
while (<GENES>) {
    chomp;
    
    # Parse out the range
    if (/^>/) {
	$_ =~ /.*:(\S*) /;
	$r = $1;
    }
    else {
	$RangeToSeq{$r} .= "$_";
    }
}


# Now have both sets of data structures.


##### 3.
# Get the average mRNA seq values for each gene and store in a hash by the ID
# Each file has 2 columns: 1. pos 2. value
# Get all the file names

my %IDtoAvg;

my $inDir = $ARGV[2];
my @files = `ls -1 $inDir/*.txt`;
foreach my $f (@files) {
    chomp $f;
    open (IN, "<$f") or die $!;
    my $total = 0;
    my $num = 0;
    while (<IN>) {
	chomp;
	my ($pos, $val) = split;
	$total += $val;
	$num++;
    }
    close IN;
    
    my $average = $total / $num;
    
    # get the ID number out of the file name string
    $f =~ /.*-(.*)-.*\.txt/;
    my $id = $1;
    
    $IDtoAvg{$id} = $average;
}


# Now have the average mRNA counts for each gene.
# Need to go through each gene that has mRNA density and find its sequence in the appropriate
# hash. Count its codons, scale them by the mRNA density, and add these to the total
# codon counts.

my %aa_table = (
             UUU => 'F', UUC => 'F', UUA => 'L', UUG => 'L', UCU => 'S',
             UCC => 'S', UCA => 'S', UCG => 'S', UAU => 'Y', UAC => 'Y',
             UAA => '*', UAG => '*', UGU => 'C', UGC => 'C', UGA => '*',
             UGG => 'W', CUU => 'L', CUC => 'L', CUA => 'L', CUG => 'L',
             CCU => 'P', CCC => 'P', CCA => 'P', CCG => 'P', CAU => 'H',
             CAC => 'H', CAA => 'Q', CAG => 'Q', CGU => 'R', CGC => 'R',
             CGA => 'R', CGG => 'R', AUU => 'I', AUC => 'I', AUA => 'I',
             AUG => 'M', ACU => 'T', ACC => 'T', ACA => 'T', ACG => 'T',
             AAU => 'N', AAC => 'N', AAA => 'K', AAG => 'K', AGU => 'S',
             AGC => 'S', AGA => 'R', AGG => 'R', GUU => 'V', GUC => 'V',
             GUA => 'V', GUG => 'V', GCU => 'A', GCC => 'A', GCA => 'A',
             GCG => 'A', GAU => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
             GGU => 'G', GGC => 'G', GGA => 'G', GGG => 'G',
             );

my @codon_list = keys %aa_table;
#All 64 codons stored in @codon_list.
 
#Hash of lists of all possible codon for each amino acid
#Input amino acid; output list of synonymous codons for that amino acid
my %all_synonyms = (
	F 	=> ['UUU', 'UUC'],
	L 	=> ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
	S 	=> ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
	Y 	=> ['UAU', 'UAC'],
        '*' 	=> ['UAA', 'UAG', 'UGA'],
	C 	=> ['UGU', 'UGC'],
	W 	=> ['UGG'],
	P       => ['CCU', 'CCC', 'CCA', 'CCG'],
	H       => ['CAU', 'CAC'],
	Q       => ['CAA', 'CAG'],
	R       => ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
	I       => ['AUU', 'AUC', 'AUA'],
	M       => ['AUG'],
	T       => ['ACU', 'ACC', 'ACA', 'ACG'],
	N       => ['AAU', 'AAC'],
	K       => ['AAA', 'AAG'],
	V       => ['GUU', 'GUC', 'GUA', 'GUG'],
	A       => ['GCU', 'GCC', 'GCA', 'GCG'],
	D       => ['GAU', 'GAC'],
	E       => ['GAA', 'GAG'],
	G       => ['GGU', 'GGC', 'GGA', 'GGG'],
	);

# Create the hash for storing the counts and initialize to zero
my %counts;

foreach (@codon_list) {
    $counts{$_} = 0;
}


foreach my $gene (keys %IDtoAvg) {
    my $avg = $IDtoAvg{$gene};
    
    my $range = $IDtoRange{$gene};
    
    my $seq = $RangeToSeq{$range};
    
    # convert sequence to RNA
    $seq =~ tr/Tt/Uu/;
    
    my @codons = $seq =~ /([AUGC]{3})/g;
    # update here: ignore first, last 10 codons,
    # since I do that for the footprint codon counts too
    for (0 .. 9) {
	pop @codons;
	shift @codons;
    }
    foreach(@codons) {
	# add AVG number of counts
	$counts{$_} += $avg;
    }
}

# Print the results suitable for rowstacked histogram (good for plotting)
# First column is name of amino acid
# Next columns are counts for synonyms of that a.a.; - for missing data

my $outFile = $ARGV[3];
open (OUT, ">$outFile") or die $!;

foreach my $aa (sort keys %all_synonyms) {
    print OUT "$aa\t";
    foreach my $c (@{$all_synonyms{$aa}}) {
	print OUT "$counts{$c}\t";
    }
    
    # Fill in the gaps with "-"
    my $numCodons = scalar @{$all_synonyms{$aa}};
    for ($numCodons .. 5) {
	print OUT "-\t";
    }
    
    print OUT "\n";
}

print OUT "\n";
close OUT;

# Run Gnuplot to generate the rowstacked plot
open GNUPLOT, "| gnuplot";
print GNUPLOT <<gnuplot_commands;
set terminal eps
set output "$outFile.eps"
set key off
set title "mRNA Codon Counts"
set ylabel "Total number in mRNA Sequencing Dataset"
set xtic nomirror
set border 11
set datafile missing "-"
set style fill solid border -1
set style data histogram
set style histogram rowstacked
plot "$outFile" u 2:xtic(1), for [i=3:7] '' using i
unset output
gnuplot_commands
close GNUPLOT;

# Print for each codon, individually
my $outFile2 = $ARGV[4];
open (OUT, ">$outFile2") or die $!;
foreach (@codon_list) {
    print OUT "$_\t$counts{$_}\n";
}
close OUT;


# Print a file that has 1. the ID and 2. the average value
my $outFile3 = $ARGV[5];
open (OUT, ">$outFile3") or die $!;
foreach (keys %IDtoAvg) {
    print OUT "$_\t$IDtoAvg{$_}\n";
}
close OUT;
