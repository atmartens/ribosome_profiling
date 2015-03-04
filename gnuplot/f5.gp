#!/usr/bin/gnuplot

set term pngcairo font "Arial, 18" size 992,645

set output "img/figure5.png"

set key off
set xlabel "Codon Position"
set ylabel "Number of Reads"
set border 3
set xtic nomirror
set ytic nomirror

set object 1 rect from 130,0 to 138,50 fc rgbcolor "red"

plot "data/amiB/amiB-b4169-density.txt" w li
unset output
