#!/usr/bin/perl -w
use strict;

# Compute the codon occurrences of ribosome footprints overall, dependent both on
# position in the footprint and categorized by the length of the footprint
# Use as input a "simplified" version of a SAM file
# Output a heatmap for each codon. Need to normalize the counts
# as percentages of total, and not absolute numbers
# 
# This version prints the heatmaps with the whitespace on the left-hand side
# (diff lengths aligned to 3' end)

# Store length of the whole genome file which is assumed to be fna format
# each genome is a member of a hash of genomes, given by the NC identifier
# Must get the identifier from the first line.
# Input a directory of genome files. Cannot contain anything else.
# 
# Drop anything > 10 codons (30 nt)

my $genomeDirectory = $ARGV[0];
my $pttDirectory = $ARGV[1];
my $fpsFile = $ARGV[2];
# directory in which to put the output files
my $out = $ARGV[3];

# max fp length, in codons
# 30 nt, for S. cerevisiae
my $maxLength = 30;
my $maxCodons = $maxLength / 3;

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
    $firstLine =~ /\|(NC_.*)\./;
    
    my $ID = $1;
    push @IDs, $ID;
    
    my $genome = "";
    while (<INGENOME>) {
        chomp;
        $genome .= $_;
    }
    
    $fnaData{$ID} = length($genome);
    
    close INGENOME;
}

# Assume that the IDs for the FNA files are the same as for the PTT files,
# use the IDs stored in @IDs
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
    for (0 .. $fnaData{$id}) {
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

# All 64 codons
my @codonlist = @{CodonList()};

# the number of occurrences
# store a hash for each length
my %counts;

# footprints file in simplified format:
open (IN, "<$fpsFile") or die $!;

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
    
    my $seq_inframe;
    
    # check every single base to make sure no frame equals zero
    foreach my $base ($startpos .. $startpos + $l) {
        if (${$pttData{$ID}}[$base] == 0) {
            # skip the whole line in the while loop
            next LINE;
        }
    }
    
    # remove 1 or 2 bases from 5' end if not multiple of 3, to keep frameness
    
    # positive strand
    if (${$pttData{$ID}}[$startpos] > 0) {
        my $remove = (4 - ${$pttData{$ID}}[$startpos]) % 3;
        $seq_inframe = substr $sequence, $remove;
    }
    
    # negative strand
    else { 
        my $remove = (4 + ${$pttData{$ID}}[$startpos + $l - 1]) % 3;
        $seq_inframe = substr $sequence, $remove;
    }
    
    my @codons = $seq_inframe =~ /([ATGC]{3})/g;
    my $length = int($l / 3);
    
    # drop anything shorter than max length
    if ($length <= $maxCodons) {
        # Check to see if we have a hash representing this length yet
        # If it doesn't, create it
        # This hash will contain n=length sub-hashes
        unless (exists($counts{$length})) {
            my @newlist;
            for my $pos (0 .. $length) {
                my %newsubhash;
                foreach my $codon (@codonlist) {
                    $newsubhash{$codon} = 0;
                }
                $newlist[$pos] = \%newsubhash;
            }
            $counts{$length} = \@newlist;
        }
        
        # Count the codons, but store in the place appropriate for FPs of this length
        # scale by the number of times that this footprint was found ~ aka in $num
        for my $pos (0 .. $#codons) {
            ${$counts{$length}}[$pos]{$codons[$pos]} += $num;
        }
    }
}

close IN;

# Need to go through the tabulations and compute, for each length,
# at each position, the total number of counts across all codons.
# Loop through the existing data structures and store these values in
# new data structures

my %sums_lengths;

foreach my $l (keys %counts) {
    my @sums;
    $sums_lengths{$l} = \@sums;
    foreach my $pos (0 .. $l - 1) {
        foreach my $c (@codonlist) {
            $sums[$pos] += $counts{$l}[$pos]{$c};
        }
    }
}


# Get the largest value in keys %counts
my @templist = sort {$a <=> $b} keys %counts;
my $largestval = $templist[$#templist];

foreach my $c (@codonlist) {
    open (OUT, ">$out/$c.heatmap.txt") or die $!;
    foreach my $l (@templist) {
      foreach my $pos (0 .. $largestval - 1) {
          if ( defined($sums_lengths{$l}[$pos]) && ($sums_lengths{$l}[$pos] > 0) ) {
              print OUT $counts{$l}[$pos]{$c} / $sums_lengths{$l}[$pos] * 100;
              print OUT "\t";
          }
          else {
              print OUT "-\t";
          }
      }
      print OUT "\n";
    }
    close OUT;
}


# Run Gnuplot a bunch of times to make plots of the data
foreach my $c (@codonlist) {
    my $file_in = "$out/$c.heatmap.txt";
    my $file_out = "$out/$c.heatmap.png";
    open GNUPLOT, "| gnuplot";
    print GNUPLOT <<gnuplot_commands;
    set terminal png
    set output "$file_out"
    set title "\\"$c\\" Codon composition in ribosome footprints"
    set ylabel "Position in Footprint"
    set xlabel "Length of Footprint"
    set autoscale
    set datafile missing "-"
    plot "$file_in" matrix with image
    unset output
gnuplot_commands
    close GNUPLOT;
}


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
