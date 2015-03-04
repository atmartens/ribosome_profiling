#!/usr/bin/perl -w
use strict;

# Given a directory containing codon position-specific tallies of mRNA sequencing read counts
# for many genes, compute the transcriptomic codon counts

# 1. read in all ptt files
# 2. extract the gene names, synonyms, and locations
# 3. read in genome FNA files (whole chromosome)
# 4. extract sequences corresponding to each PTT entry's location
# 5. go through each pre-computed 1D mRNAseq read density file (density / codon) and compute
#    an average density for each gene == mRNA level
# 6. count codons in each sequence and scale it by its mRNA level, track the sum
# 7. print codon counts to file -> these are the codon counts weighted by mRNA expression

# ARGUMENTS:
# 0. PTT files directory
# 1. FNA files directory
# 2. mRNA seq directory
# 3. output file for each codon alone
# 4. output file for rowstacked histogram
# # 5. output file with ID and avg value

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

# Read in the PTT files

my $pttDirectory = $ARGV[0];
my $genomeDirectory = $ARGV[1];
my $RNAseqDirectory = $ARGV[2];
my $countsOutputFile = $ARGV[3];
my $histoOutFile = $ARGV[4];


my @pttFiles = `ls -1 $pttDirectory`;

my %IDtoPTT;

foreach my $pttFile (@pttFiles) {
    chomp $pttFile;
    
    my $inFile = "$pttDirectory/$pttFile";
    
    # assume file ends in .ptt, remove this to get the ID
    my $l = length($pttFile);
    my $ID = substr $pttFile, 0, $l - 4;
    
    # open file, drop 3 lines, and dump contents into hash
    open (PTT, "<$inFile") or die $!;
    
    # Drop the first 3 lines
    for (0 .. 2) {
        <PTT>;
    }
    
    # Keep the rest
    my @pttLines;
    
    while (<PTT>) {
	chomp;
	push @pttLines, $_;
    }
    
    close PTT;
    
    $IDtoPTT{$ID} = \@pttLines;
}

my @fnaFiles = `ls -1 $genomeDirectory`;

# Hash to store each genome FNA file
my %fnaData;

# Get the genome info, store in hash
foreach my $fna (@fnaFiles) {
    chomp $fna;
    my $inFile = "$genomeDirectory/$fna";
    open (INGENOME, "<$inFile") or die $!;    
    
    # get the identifier
    my $firstLine = <INGENOME>;
    # sometimes can get issues with this Regex: diff files, diff patterns
    $firstLine =~ /\|(NC_[^.]*)/;
    my $ID = $1;
    
    # get whole chromosome seq.
    my $genome = "";
    while (<INGENOME>) {
        chomp;
        $genome .= $_;
    }
    
    $fnaData{$ID} = $genome;
    
    close INGENOME;
}


# iterate through PTT, do genomic lookups to get sequences and do codon counts
# simultaneously look up and compute mRNA read density


# Store codon counts for all genes, weighted by expression level
# initialize to zero
my %codonCounts;
foreach (keys %aa_table) {
    $codonCounts{$_} = 0;
}

# something wrong with minus strand...

foreach my $chromosomeID (keys %IDtoPTT) {
    foreach my $pttLine (@{$IDtoPTT{$chromosomeID}}) {
	# get location
	my @line = split /\t/, $pttLine;
	
	my $range = $line[0];
	my ($left, $right) = split /\.\./, $range;
	
	my $strand = $line[1];
	
	my $gene = $line[4];
	my $synonym = $line[5];
	
	# Look up range in FNA
	my $seq = substr $fnaData{$chromosomeID}, $left - 1, $right - $left + 1;
	
	# if minus strand, make revcomp
	if ($strand eq "-") {
	    $seq =~ tr/ATGCatgc/TACGtacg/;
	    $seq = reverse $seq;
	}
	
	# convert to RNA
	$seq =~ tr/Tt/Uu/;
	
	# chop up into codons
	my @codons = $seq =~ /([AUGC]{3})/g;
	
	###
	
	# Input file with read density, formatted in 2 columns: 1. codon pos 2. num
	# compute average
	my $infile = "$RNAseqDirectory/$gene-$synonym-density.txt";
	
	# check that it exists (some genes are not expressed):
	if (-e $infile) {
	    open (IN, "<$infile") or die $!;
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
	    
	    ###
	    
	    # tally:
	    foreach my $c (@codons) {
		$codonCounts{$c} += $average;
	    }
	}
    }
}
 
# Hash of lists of all possible codon for each amino acid
# Input amino acid; output list of synonymous codons for that amino acid
# Useful for rowstacked hisotgram
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


# Print the results suitable for rowstacked histogram (good for plotting)
# First column is name of amino acid
# Next columns are counts for synonyms of that a.a.; - for missing data

open (OUT, ">$histoOutFile.tsv") or die $!;

foreach my $aa (sort keys %all_synonyms) {
    print OUT "$aa\t";
    foreach my $c (@{$all_synonyms{$aa}}) {
	print OUT "$codonCounts{$c}\t";
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
set output "$histoOutFile.eps"
set key off
set title "mRNA Codon Counts"
set ylabel "Total number in mRNA Sequencing Dataset"
set xtic nomirror
set border 11
set datafile missing "-"
set style fill solid border -1
set style data histogram
set style histogram rowstacked
plot "$histoOutFile.tsv" u 2:xtic(1), for [i=3:7] '' using i
unset output
gnuplot_commands
close GNUPLOT;

# Print for each codon, individually
open (OUT, ">$countsOutputFile") or die $!;
foreach (keys %aa_table) {
    print OUT "$_\t$codonCounts{$_}\n";
}
close OUT;


# # # Print a file that has 1. the ID and 2. the average value
# # my $outFile3 = $ARGV[5];
# # open (OUT, ">$outFile3") or die $!;
# # foreach (keys %IDtoAvg) {
# #     print OUT "$_\t$IDtoAvg{$_}\n";
# # }
# # close OUT;


