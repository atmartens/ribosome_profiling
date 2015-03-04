#!/usr/bin/gnuplot

set term pngcairo font "Arial, 18" size 992,645

set datafile missing "-"
set xrange [-0.5:46.5]
set yrange [20.5:47.5]
set key off
set cblabel "GC Content"
set xlabel "Position in read (nt)"
set ylabel "Read length (nt)"

set output "img/figure9.png"
plot "data/rnaseq-gc/gc_LtoR.csv" matrix with image
unset output
