#!/usr/bin/perl -w
use strict;
use PDL;
use PDL::ImageND;

# Take a DNA or RNA sequence and compute a window-averaged GC content
# using the Savitzky-Golay smoothing algorithm

# fasta; assume the seq is just on one line, right after the > line
my $sequenceFile = $ARGV[0];

open (IN, "<$sequenceFile") or die $!;
<IN>;
my $seq = <IN>;
close IN;

# cut up sequence into nucleotides in a list
my @list_nuc = $seq =~ /[AUTGC]/g;

# convert to GC: 1 if G or C, else 0
my @list_one_or_zero;
foreach (@list_nuc) {
    if (/[AUT]/) {
        push @list_one_or_zero, 0;
    }
    else {
        push @list_one_or_zero, 1;
    }
}

# SavGol params
my $nleft = 13;
my $nright = 13;
my $degree = 4;
my $SGFilter = SavGol($nleft,$nright,$degree);

my @SGData = runSavGol($SGFilter, @list_one_or_zero);

for (my $i = 1; $i < scalar(@SGData) + 1; $i++) {
    print "$i\t$SGData[$i-1]\n";
}



# Precompute a Savitzy Golay matrix
# Arguments:    1. numleft (nl)         where window size = numleft + numright + 1
#               2. numright (nr)
#               3. degree of the polynomial (deg) (usually 4)
#               Returns a piddle with which to filter data.
sub SavGol {
    my $nl = shift;
    my $nr = shift;
    my $deg = shift;
    
    my $i = $nl+$nr+1;
    
    # Start making the matrix
    my $x = ((PDL->zeroes($i)->xvals) - $nl)->float;
    my $AT = ((PDL->zeroes($i,$deg+1)->xvals))->float;
    
    for(0..$deg) {
        (my $tmp = $AT->slice(":,($_)")) .= ($x ** $_);
    }
    
    # (AT x A)^-1 x AT
    # grab the top row
    return transpose(transpose(inv($AT x transpose($AT)) x $AT)->slice("0:0"));
}


# Given a Savitzy Golay coefficient matrix, and some data, convolve the two
# Arguments:    1. SavGol filter
#               2. Data
sub runSavGol {
    my $SGFilter = shift;
    
    # convert to piddle & make it a matrix with the proper dimensions
    # run the convolution
    # convert back to a regular perl array & return
    return list convolveND(pdl(@_)->dummy(1), $SGFilter);
}
