#!/usr/bin/perl -w
use strict;

# Compute the codon occurrences of ribosome footprints overall, dependent both on
# position in the footprint and categorized by the length of the footprint
# Use as input a "simplified" version of a SAM file
# Output a heatmap for each codon. Normalize the counts
# as percentages of total, and not absolute numbers
# 
# Store length of the whole genome file which is assumed to be fna format
# each genome is a member of a hash of genomes, given by the NC identifier
# Must get the identifier from the first line.
# Input a directory of genome files. Cannot contain anything else.
# 
# fill in seqs on 5' end, align to 3' end
# 

my $genomeDirectory = $ARGV[0];
my $pttDirectory = $ARGV[1];
my $fpsFile = $ARGV[2];
my $outDir = $ARGV[3];
my $lengthTalliesFile = $ARGV[4];

# max fp length, in codons
my $maxLength = 45;
my $maxCodons = int($maxLength / 3);

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

# All 64 codons
my @codonlist = @{CodonList()};

# the number of occurrences
# store a hash for each length
my %counts;

# footprints file in simplified format:
open (IN, "<$fpsFile") or die $!;

# Store the number in each length category
my %lenCat;

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
    
    # keep multiple of 3 for frame at the end:
    # remove 1 or 2 nucleotides at the very end
    while ( ($l % 3) != 0 ) {
        $l--;
    }
    
    # Look up the actual sequence from the genome and fill in the gaps
    # need to preserve frame information: adjust accordingly by shifting towards 5'
    
    # positive strand
    my $fullSeq;
    if (${$pttData{$ID}}[$startpos] > 0) {
        while ( ${$pttData{$ID}}[$startpos] != 1 ) {
            $startpos--;
        }
        $fullSeq = substr $fnaData{$ID}, ($startpos - 1) - ($maxLength - $l), $maxLength;
    }
    
    # negative strand
    else {
        while (${$pttData{$ID}}[$startpos] != -3 ) {
            $startpos++;
        }
        $fullSeq = substr $fnaData{$ID}, ($startpos - 1), $maxLength;
        $fullSeq =~ tr/ATGC/TACG/;
        $fullSeq = reverse $fullSeq;
    }
    
    my @codons = $fullSeq =~ /([ATGC]{3})/g;
    
    my $length = int($l / 3);
    
    # Tally the length categories ~ total number of FPs of each length (in codons)
    $lenCat{$length} += $num;
    
    # Check to see if we have a hash representing this length yet
    # If it doesn't, create it
    # This hash will contain n=15 sub-hashes (45 nt --> 15 codons)
    unless (exists($counts{$length})) {
        my @newlist;
        for my $pos (0 .. ($maxLength / 3) - 1 ) {
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

close IN;

# Need to go through the tabulations and compute, for each length,
# at each position, the total number of counts across all codons.
# Loop through the existing data structures and store these values in
# new data structures

my %sums_lengths;

foreach my $l (keys %counts) {
    my @sums;
    $sums_lengths{$l} = \@sums;
    foreach my $pos (0 .. $maxCodons) {
	foreach my $c (@codonlist) {
	    $sums[$pos] += $counts{$l}[$pos]{$c};
	}
    }
}

# Print a file containing the number in each length category
# and ditch anything less than some minimum number
# because of low sample size noise
open (OUT, ">$lengthTalliesFile") or die $!;
foreach my $i (sort {$a <=> $b} keys %lenCat) {
    print OUT "$i\t$lenCat{$i}\n";
    # if under 30,000 total FPs of this length category
    # then drop from the %counts hash
    if ($lenCat{$i} < 30000)  {
        delete $counts{$i};
    }
}
close OUT;

# Get the largest value in keys %counts
my @templist = sort {$a <=> $b} keys %counts;

# For Gnuplot Y ranges
my $smallest_Y_Val = $templist[0] - 0.5;
my $largest_Y_Val = $templist[$#templist] + 0.5;

# For Gnuplot X ranges
my $smallest_X_Val = -0.5;
my $largest_X_Val = $maxCodons - 0.5;

# Print the results in a format good for gnuplot heatmaps
foreach my $c (@codonlist) {
    open (OUT, ">$outDir/$c.heatmap.txt") or die $!;
    
    # First, spit out a bunch of dashes for empty lines
    for (0 .. $templist[0] - 1) {
        print OUT "-";
        for (0 .. $maxCodons - 2) {
            print OUT "\t-";
        }
        print OUT "\n";
    }
    
    foreach my $l (@templist) {
        # Need to re-do the order for the x-plotting
        # so that the right-most value is the 5' end,
        # but then we plot the x axis reversed so that we
        # see the right-most as the 3' end
        for (my $pos = $maxCodons - 1; $pos >= 0; $pos--) {
            if ( defined($sums_lengths{$l}[$pos]) && ($sums_lengths{$l}[$pos] > 0) ) {
                print OUT $counts{$l}[$pos]{$c} / $sums_lengths{$l}[$pos] * 100;
            }
            else {
                print OUT "-";
            }
            print OUT "\t";
        }
        print OUT "\n";
    }
    close OUT;
}

# Run Gnuplot a bunch of times to make plots of the data
foreach my $c (@codonlist) {
    my $U_version = $c;
    $U_version =~ tr/Tt/Uu/;
    my $file_in = "$outDir/$c.heatmap.txt";
    my $file_out = "$outDir/$c.heatmap.png";
    open GNUPLOT, "| gnuplot";
    print GNUPLOT <<gnuplot_commands;
    set terminal png
    set key off
    set output "$file_out"
    set ylabel "Length of footprint (codons)"
    set xlabel "Position in footprint"
    set cblabel "Codon $U_version %"
    set autoscale
    set datafile missing "-"
    set yrange [$smallest_Y_Val:$largest_Y_Val]
    set xrange [$largest_X_Val:$smallest_X_Val]
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
