#!/usr/bin/gnuplot
set terminal pngcairo font "Arial, 18" size 992,645
set key off
set output "img/figure6.png"
set multiplot layout 2,2

set linestyle 1 lt 3 lw 4 linecolor rgb "red"

set label 2 "5\'" at graph -0.05,-0.055 font ',14'
set label 3 "3\'" at graph 1,-0.055 font ',14'

set xrange [-0.5:10.5]
set yrange [6.5:12.5]
set datafile missing "-"


set cbrange [0:]
set cblabel offset 1.5

set xtics scale 0
set xtics 0, 2, 10
set xtics offset 0,-0.3
set ytics scale 0

set ytics format "%.0f"
set ytics 0, 1, 13

set title "CAC / 3-AT"
set cblabel "% of all codons"
set cbtics 1 format "%.0f"
set ylabel "Footprint length (codons)" offset 0,-6
plot "data/yeast_his/CAC.3at.csv" matrix with image

set title "CAU / 3-AT"
set cblabel "% of all codons"
set cbtics 1 format "%.0f"
set ylabel " "
set xlabel " "
plot "data/yeast_his/CAT.3at.csv" matrix with image

set yrange [6.5:11.5]

set title "CAC / untreated"
set cblabel "% of all codons"
set cbtics 1 format "%.0f"

set xlabel "Distance from 5' end (codons)" offset 20,0
plot "data/yeast_his/CAC.untreated.csv" matrix with image

set title "CAU / untreated"
set cblabel "% of all codons"
set cbtics 1 format "%.0f"
set xlabel " "
plot "data/yeast_his/CAT.untreated.csv" matrix with image

unset multiplot
