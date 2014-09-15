#!/usr/bin/perl -w
use strict;

# Go through a PTT file and remove any elements which have apparent
# internal stop codons. Print a new, reduced, PTT file.

my %aa_table = (
             TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L', TCT => 'S',
             TCC => 'S', TCA => 'S', TCG => 'S', TAT => 'Y', TAC => 'Y',
             TAA => '*', TAG => '*', TGT => 'C', TGC => 'C', TGA => '*',
             TGG => 'W', CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',
             CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P', CAT => 'H',
             CAC => 'H', CAA => 'Q', CAG => 'Q', CGT => 'R', CGC => 'R',
             CGA => 'R', CGG => 'R', ATT => 'I', ATC => 'I', ATA => 'I',
             ATG => 'M', ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',
             AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K', AGT => 'S',
             AGC => 'S', AGA => 'R', AGG => 'R', GTT => 'V', GTC => 'V',
             GTA => 'V', GTG => 'V', GCT => 'A', GCC => 'A', GCA => 'A',
             GCG => 'A', GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
             GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',
             );

# Store the whole genome file which is assumed to be fna format
# each genome is a member of a hash of genomes, given by the NC identifier
# Must get the identifier from the first line.
# Input a directory of genome files. Cannot contain anything else.

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


my $pttDirectory = $ARGV[1];
my $new_pttDirectory = $ARGV[2];

my @pttFiles = `ls -1 $pttDirectory`;

# Go through PTT files and remove anything with internal stop codons
foreach my $ptt (@pttFiles) {
    chomp $ptt;
    my $inFile = "$pttDirectory/$ptt";
    my $outFile = "$new_pttDirectory/$ptt";
    
    # assume file ends in .ptt, remove this to get the ID
    my $l = length($ptt);
    my $ID = substr $ptt, 0, $l - 4;
    
    open (PTT, "<$inFile") or die $!;
    open (NEWPTT, ">$outFile") or die $!;
    
    # Keep the first 3 lines
    for (0 .. 2) {
        my $temp = <PTT>;
        print NEWPTT $temp;
    }
    
    while (<PTT>) {
        chomp;
        my @line = split /\t/;
        
        my $range = $line[0];
        my $strand = $line[1];
        
        # parse the range to get the beginning and end
        my ($left, $right) = split /\.\./, $range;
        
        my $seq = substr $fnaData{$ID}, $left - 1, $right - ($left - 1);
        if ($strand eq "-") {
            $seq = reverse $seq;
            $seq =~ tr/ATGC/TACG/;
        }
        
        # Translate, remove final stop codon, see if there are any others
        # If not, go ahead and keep this PTT line
        # convert to codons
        my @codons = $seq =~ /([ATGC]{3})/g;
        
        # remove last codon, presumably a stop
        pop @codons;
        
        # translate
        my @aa;
        foreach my $c (@codons) {
            push @aa, $aa_table{$c};
        }
        
        # count number of matches to *
        my $num_times = grep /\*/, @aa;
        
        if ($num_times == 0) {
            print NEWPTT "$_\n";
        }
    }
    
    close PTT;
    close NEWPTT;
}
