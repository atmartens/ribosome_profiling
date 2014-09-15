#!/usr/bin/perl -w
use strict;

# Go through a SAM file and "fix" the sequences so that they are perfect matches to the reference
# i.e. swap all Ns to the reference base, and all other mismatches
# output a simplified file with only the sequences, their mapped positions and their occurrences
# 
# This version keeps chromsome information. i.e. which NC_ file the info came from

# Store the whole genome file which is assumed to be fna format
# each genome is a member of a hash of genomes, given by the NC identifier
# Must get the identifier from the first line.
# Input a directory of genome files. Cannot contain anything else.
# 
# Drop anything > 50 nt, assuming that's the max read length.

my $maxLength = 50;

my $genomeDirectory = $ARGV[0];

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

# footprints file in SAM format
# assume that the comments with @ are excluded already
open (IN, "<$ARGV[1]") or die $!;

my %counts;

while (<IN>) {
    # $line[0] : the ID of the read
    # $line[1] : the bit score ~ 4 means did not align, 16 is complementary strand
    # $line[2] : the ID of the genome, including the NC file name
    # $line[3] : the chromosomal position of the alignment
    # $line[9] : the sequence of the read
    chomp;
    my @line = split /\t/;
    
    my $org = $line[2];
    
    $org =~ /\|(NC_.*)\./;
    my $ID = $1;
    
    my $sequence = $line[9];
    my $startpos = $line[3];
    my $length = length($sequence);
    
    if ($length <= $maxLength) {
        my $reference = substr $fnaData{$ID}, $startpos - 1, $length;
        
        if ($line[1] == 16) {
            $reference =~ tr/ATGCatgc/TACGtacg/;
            $reference = reverse $reference;
        }
        
        $reference = "$reference\t$startpos\t$ID";
        
        $counts{$reference}++;
    }
}


close IN;

open (OUT, ">$ARGV[2]") or die $!;

foreach (keys %counts) {
    print OUT "$_\t$counts{$_}\n";
}

close OUT;
