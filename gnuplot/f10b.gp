#!/usr/bin/gnuplot

# set term pdf size 8,6
set term pngcairo font "Arial, 18" size 992,645

set output "img/f10b.png"

set datafile missing "-"

set key off
set xlabel "3' end of footprint (nucleotides)"
set ylabel "5' end of footprint (nucleotides)"

set cblabel "Number of footprints"
set cblabel offset character 1,0

set mxtics 4
set mytics 4

set xtics out
set ytics out

set xrange [367.5:472.5]
set yrange [338.5:442.5]

# so that it doesn't show black as zero... cause then what's white?
set cbrange [0.1:]

plot "data/amiB/amiB-b4169.heat.txt" matrix w image
unset output
