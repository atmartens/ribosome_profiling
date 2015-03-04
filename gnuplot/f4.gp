#!/usr/bin/gnuplot
set terminal pngcairo font "Arial, 18" size 992,645
set key off
set output "img/figure4.png"
set multiplot layout 2,2 rowsfirst

set linestyle 1 lt 3 lw 4 linecolor rgb "red"

set macros
POS = "at graph -0.05,1.02 font ',20'"

set label 2 "5\'" at graph -0.05,-0.05 font ',14'
set label 3 "3\'" at graph 1,-0.05 font ',14'

set xrange [9.5:-0.5]
set yrange [6.5:10.5]
set datafile missing "-"


unset colorbox
set cbrange [0:2.5]

set xtics scale 0
set ytics scale 0

set ytics format "%.0f"
set ytics 0, 1, 10

set xlabel " "

set title "CCA"
set ylabel "Footprint length (codons)" offset 0,-5
plot "data/ecoli_pro/CCA.heatmap.txt" matrix with image

set ylabel " "

set title "CCU"
plot "data/ecoli_pro/CCT.heatmap.txt" matrix with image


set title "CCC"
plot "data/ecoli_pro/CCC.heatmap.txt" matrix with image

set title "CCG"

set colorbox horizontal user origin 0.385,0.475 size .25,.025
set cblabel "Fold Enrichment" offset 0,4.15
set cbtics offset 0,0.3

set termoption enhanced
# set label 4 "E. coli proline codon enrichment" at 15.85,20.5

set xlabel "Distance from 3' end (codons)" offset -18,0

plot "data/ecoli_pro/CCG.heatmap.txt" matrix with image

unset multiplot
