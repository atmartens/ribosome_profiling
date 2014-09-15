#!/usr/bin/perl -w
use strict;

# Go through a "simplified" SAM file and discard footprints which are not in protein coding
# regions or which are within the first or last 10 codons of protein coding regions
# Output a file of the same format
# By storing things in memory it use several gigabytes, but noticeable speedup with sorting
# # (roughly 1 hour versus 17 hours)


my $pttDirectory = $ARGV[0];
my $samFile = $ARGV[1];
my $outFile = $ARGV[2];

my @pttFiles = `ls -1 $pttDirectory`;

# input ID, output a list of PTT lines
my %ptt;

foreach my $nextFile (@pttFiles) {
    chomp $nextFile;
    
    # assume ID = the ptt filename before extension ".ptt"
    my $l = length($nextFile);
    my $id = substr $nextFile, 0, $l - 4;
    
    my @newlist;
    
    open (IN, "<$pttDirectory/$nextFile") or die $!;
    
    for (0 .. 2) {
        <IN>;
    }
    
    while (<IN>) {
        chomp;
        push @newlist, $_;
    }
    
    close IN;
    
    $ptt{$id} = \@newlist;
}

# Read in the simplified SAM file: 1. seq 2. pos 3. chromo id 4. num
open (SAM, "<$samFile") or die $!;
my @simple_SAM = <SAM>;
close SAM;

foreach (@simple_SAM) {
    chomp;
}

# sort it by the 2nd column i.e. the mapped position
# using the "map sort map" technique (Schwartzian Transform)
@simple_SAM = map {$_->[0]}
             sort { $a->[2] <=> $b->[2] }
             map {chomp;[$_,split(/\t/)]} @simple_SAM;

# Go through the sorted list and compare it to the PTT elements
# until we find a PTT element which contains the position
# and the position is further than 10 codons from the extremities

open (OUT, ">$outFile") or die $!;

foreach my $sam_info(@simple_SAM) {
    chomp $sam_info;
    my ($seq, $pos, $id, $num) = split /\t/, $sam_info;
    
    my $front = $pos;
    my $back = $pos + length($seq);
    
    # @{$ptt{$id}} is the de-referenced array reference to the PTT list corresponding
    # to the chromosome $id stored in hash %ptt
    
    my $num_ptt = scalar(@{$ptt{$id}});
    PTT: for (my $i = 0; $i < $num_ptt; $i++) {
	# location, strand, length, PID, gene, synonym, code, COG, product
	my @line = split /\t/, @{$ptt{$id}}[$i];
	my $range = $line[0];
	my ($left, $right) = split /\.\./, $range;
	
	# Because SAM file is sorted, if this PTT element has a smaller position, then
	# we can remove it from the list, because nothing will ever cross it again.
	if ( ($left < $front) && ($right < $front) && ($left < $back) && ($right < $back) ) {
	    # remove this element
	    splice @{$ptt{$id}}, $i, 1;
	    # minus one, or else we go out of bounds
	    $num_ptt--;
	    # minus one, since we just eliminated an element
	    $i--;
	}
	
	# The footprint is within the bounds of this gene:
	# keep it and skip all the other genes, going to the next footprint
	if (($front >= $left + 30 ) && ($back <= $right - 30)) {
	    print OUT "$sam_info\n";
	    last PTT;
	}
    }
}

close OUT;
