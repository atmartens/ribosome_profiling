#!/usr/bin/perl -w
use strict;

# Given a file (or more), for a gene, with the beginning, endings of the
# footprints:
# Produce a file containing a 2D matrix of densities of the footprints
# need to get the length of the gene, based on its name, from the ptt file
# 
# file containing gene coordinates + names: use a PTT file

open (COORDINATES, "<$ARGV[0]") or die $!;
my @coords = <COORDINATES>;
close COORDINATES;

# drop the first 3 lines of the PTT file
for (0 .. 2) {
    shift @coords;
}

my $inDir = $ARGV[1];
my $outDir = $ARGV[2];

foreach (@coords) {
    chomp;
    
    my @line = split /\t/;
    # location, strand, length, PID, gene, synonym, code, COG, product
    
    my $range = $line[0];
    
    # parse the range to get the beginning and end
    my ($min, $max) = split /\.\./, $range;
    
    my $strand = $line[1];
    my $length = $line[2] * 3;
    my $name = $line[4];
    my $synonym = $line[5];
    my $prod = $line[8];
    
    my $name_synonym = "$name-$synonym";
    
    # Create the list of hashes
    my @dim1;
    for my $i (0 .. $length) {
	my %dim2;
	$dim1[$i] = \%dim2;
    }
    
    # These files contain the beginning and ending coordinates of the reads for a gene
    # and have already been generated (part 1) and then combined
    my $infile = "$inDir/$name_synonym.txt";
    
    # Some of these won't exist because no reads were mapped...
    if (-e $infile) {
        open (IN, "<$infile") or die $!;
        my $counter = 0;
        while (<IN>) {
            chomp;
            my ($left, $right, $num) = split;
            
            # Data are in nucleotides
            # Need to convert from chromosome positions to
            # relative pos. from start codon
            
            # make sure that + or - strand genes are treated correctly
            my $left_abs;
            my $right_abs;
            
            if ($strand eq "+") {
                  $left_abs = $left - $min;
                  $right_abs = $right - $min;
            }
            else {
                  $left_abs = -1 * ($right - $max);
                  $right_abs = -1 * ($left - $max);
            }
            
            # array positions represent 5' ends
            # array values represent hashes
            # the hashes' keys represent 3' ends
            # the values are the number of footprints at the combined
            # 5'-3' positions
            $dim1[$left_abs]{$right_abs} = $num;
        }
        close IN;

	my $outfile = "$outDir/$name_synonym.heat.txt";
	open (OUT, ">$outfile") or die $!;
	# Go through each hash in the list
	for my $i (0 .. $length) {
	   
	   # get all the "ends" from the hash, sorted by order, ascending
	   my @ends = sort {$a <=> $b} keys %{$dim1[$i]};
	   
	   # make sure the hash isn't empty ~ that it has keys
	   if (scalar(@ends) > 0) {
		
		# print blank chars up to the first (lowest) key of the hash
		for my $j (0 .. $ends[0] - 1) {
		    print OUT "-\t";
		}
		
		# print the values
		for my $j ($ends[0] .. $ends[$#ends]) {
		    # there can still be gaps, so make sure that each one exists
		    if (exists( ${$dim1[$i]}{$j}) ) {
			print OUT "${$dim1[$i]}{$j}";
		    }
		   
		    # this element wasn't present, so make a blank char
		    else {
			print OUT "-";
		    }
		    
		    print OUT "\t"
		}
		
		# fill in the remaining blanks
		for my $j ($ends[$#ends]+1 .. $length) {
		    print OUT "-\t";
		}
	    }
	    
	    # the hash is empty, so just print a bunch of blank chars, aka "-"
	    else {
		for my $j (0 .. $length) {
		    print OUT "-\t";
		}
	    }
	    
	    print OUT "\n";
	}
	close OUT;
    }
}
