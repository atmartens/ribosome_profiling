#!/usr/bin/gnuplot

set terminal pngcairo font "Arial, 18" size 992,645 dashed

set linestyle 1 lt 3 lw 2 linecolor rgb "white"
set key off
set samples 10000

set output "img/figure7.png"
set xlabel "Distance from 3' end (nucleotides)"
set ylabel "Footprint length (nucleotides)"
set cblabel "GC Content"
set xrange [45.5:0.5]
set yrange [20.5:45.5]
set xtic nomirror
set ytic nomirror
plot "data/gc/ecoli.csv" matrix with image, (21.49 <= x && x < 45.45) ? floor(x+0.5)-0.5 : 0/0 ls 1
unset output 
