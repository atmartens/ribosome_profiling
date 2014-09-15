#!/usr/bin/perl -w
use strict;

# Take a simplified SAM file which has been extended 5' (to 45 nt or so)
# Group the footprints by length
# Count the ATGC by position
# Make some plots

open (IN, "<$ARGV[0]") or die $!;

my %lengths;
my @letters = ("A","T","G","C");
my $extended_length;

while (<IN>) {
    chomp;
    my @line = split;
    # 1. seq
    my $s = $line[0];
    # 2. chromosome ID
    my $id = $line[1];
    # 3. number
    my $n = $line[2];
    # 4. original length
    my $l = $line[3];
    
    $extended_length = length($s);
    
    # First time reaching a footprint with this original length:
    # need to initialize data structures
    unless (defined($lengths{$l})) {
	my @positions;
	for my $counter (0 .. $extended_length - 1) {
	    my %chars;
	    for my $letter (@letters) {
		$chars{$letter} = 0;
	    }
	    $positions[$counter] = \%chars;
	}
	$lengths{$l} = \@positions;
    }
    
    my @array = $s =~ /([ATGC]{1})/g;
    
    for my $counter (0 .. $extended_length - 1) {
	$lengths{$l}[$counter]{$array[$counter]} += $n;
    }
}

close IN;

# Done tallying through the file
# Now can loop through the data structure and make some plots

foreach my $l (keys %lengths) {
    # First, need to write each length dataset to its own file
    # Compute percentages
    my $total = 0;
    
    # assume total num is same at each position
    foreach my $letter (@letters) {
	$total += $lengths{$l}[0]{$letter};
    }
    
    my $file_out = "$ARGV[1]/result-$l";
    
    open (OUT, ">$file_out.tsv") or die $!;
    print OUT "Position\tA\tU\tG\tC\n";
    for my $n (0 .. $extended_length - 1) {
        print OUT "$n\t";
        print OUT "$lengths{$l}[$n]{'A'}\t";
        print OUT "$lengths{$l}[$n]{'T'}\t";
        print OUT "$lengths{$l}[$n]{'G'}\t";
        print OUT "$lengths{$l}[$n]{'C'}\n";
    }
    close OUT;
    
    open GNUPLOT, "| gnuplot";
    print GNUPLOT <<gnuplot_commands;
    set terminal png
    set output "$file_out.png"
    set key bottom left
    set title "Length $l"
    set xlabel "Position in 5'-extended footprint"
    set ylabel "Nucleotide Frequency"
    set border 11
    set xtic nomirror
    plot for [i = 2 : 5] "$file_out.tsv" u 1:i ti col w li
    unset output
gnuplot_commands
    close GNUPLOT;
}
