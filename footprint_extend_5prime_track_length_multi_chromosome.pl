#!/usr/bin/perl -w
use strict;

# Take a "simplified" set of footprints and "complete" them such that they are all equal to 45 nt.
# If they were longer than 45 nt, drop them (very tiny fraction)
# Drop footprints which have frame of zero ~ aka which overlap non-coding regions
# 
# Keep track of the original length
# This version works with multiple chromosomes

# Store length of the whole genome file which is assumed to be fna format
# each genome is a member of a hash of genomes, given by the NC identifier
# Must get the identifier from the first line.
# Input a directory of genome files. Cannot contain anything else.

my $genomeDirectory = $ARGV[0];

my $maxLength = 45;

my @fnaFiles = `ls -1 $genomeDirectory`;
my @IDs;

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
    push @IDs, $ID;
    
    my $genome = "";
    while (<INGENOME>) {
	chomp;
	$genome .= $_;
    }
    
    $fnaData{$ID} = $genome;
    
    close INGENOME;
}

# Assume that the IDs for the FNA files are the same as for the PTT files,
# use the IDs stored in @IDs

my $pttDirectory = $ARGV[1];

# hash for each chromosome: each element will store a list; each element of the list is the
# frame info for that position
my %pttData;

foreach my $id (@IDs) {
    my $fileName = "$pttDirectory/$id.ptt";
    
    open (PTT, "<$fileName") or die $!;
    # coordinate info for the current PTT file
    my @coords = <PTT>;
    close PTT;

    # drop the first 3 lines of the PTT file
    for (0 .. 2) {
	shift @coords;
    }

    # store the frame information; initialize to zero
    my @frames;
    for (0 .. length($fnaData{$id})) {
	$frames[$_] = 0;
    }

    # annotate across the length of the chromosome,
    # for each position, whether it is in the 1st, 2nd, or 3rd position
    # in the codon (aka the frame)
    # the entire genome of the organism
    foreach (@coords) {
	# location, strand, length, PID, gene, synonym, code, COG, product
	my @line = split /\t/;
	
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
    
    # Store the @frames into the hash
    $pttData{$id} = \@frames;
}

# footprints file in simplified format
open (IN, "<$ARGV[2]") or die $!;

# write to this file the "extended" sequences
open (OUT, ">$ARGV[3]") or die $!;

LINE: while (<IN>) {
    # sequence, mapped position, counts
    # need to determine if frame is correct
    chomp;
    my @line = split /\t/;
    
    my $sequence = $line[0];
    my $startpos = $line[1];
    my $ID = $line[2];
    my $num = $line[3];
    my $l = length($sequence);
    
    # Skip anything > 45 nt.
    if ($l > $maxLength) {
	next LINE;
    }
    
    my $seq_inframe;
    
    # check every single base to make sure no frame equals zero
    foreach my $base ($startpos .. $startpos + $l) {
	if (${$pttData{$ID}}[$base] == 0) {
	    # skip the whole line in the while loop
	    next LINE;
	}
    }
    
    # Look up the actual sequence from the genome and fill in the gaps
    
    # positive strand
    my $fullSeq;
    if (${$pttData{$ID}}[$startpos] > 0) {
	$fullSeq = substr $fnaData{$ID}, ($startpos - 1) - ($maxLength - $l), $maxLength;
    }
    
    # negative strand
    else {
	$fullSeq = substr $fnaData{$ID}, ($startpos - 1), $maxLength;
	$fullSeq =~ tr/ATGC/TACG/;
	$fullSeq = reverse $fullSeq;
    }
    
    # Print the seq., the ID, the number of times and the original length
    print OUT "$fullSeq\t$ID\t$num\t$l\n";
    
}

close IN;
close OUT;
