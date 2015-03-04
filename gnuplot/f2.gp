#!/usr/bin/gnuplot
set terminal pngcairo font "Arial, 18" size 2102,976 dashed
set samples 10000
set key off
set output "img/figure2.png"
set multiplot layout 2,3

set linestyle 1 lt 3 lw 4 linecolor rgb "red"

set label 2 "5\'" at graph -0.05,-0.03 font ',14'
set label 3 "3\'" at graph 1,-0.03 font ',14'

set xrange [14.5:-0.5]
set yrange [7.5:13.5]
set datafile missing "-"

set cbrange [0:]
set cblabel offset 1.5

set xtics scale 0
set xtics 0, 3, 15
set xtics offset 0,-0.3

set ytics scale 0
set ytics format "%.0f"
set ytics 0, 1, 13

set xlabel " "

set cblabel "% UCA"
set cbtics 0.5 format "%.1f"
set ylabel "Footprint length (codons)" offset 0, -9.5
plot "data/ecoli_ser/TCA.heatmap.txt" matrix with image, (7.49 <= x && x < 12.55) ? floor(x+1.5)-0.5 : 0/0 ls 1

set cblabel "% UCU"
set cbtics 1 format "%.0f"
set ylabel " "
plot "data/ecoli_ser/TCT.heatmap.txt" matrix with image, (7.49 <= x && x < 12.55) ? floor(x+1.5)-0.5 : 0/0 ls 1

set cblabel "% UCC"
set cbtics 1 format "%.0f"
plot "data/ecoli_ser/TCC.heatmap.txt" matrix with image, (7.49 <= x && x < 12.55) ? floor(x+1.5)-0.5 : 0/0 ls 1

set cblabel "% UCG"
set cbtics 0.5 format "%.1f"
plot "data/ecoli_ser/TCG.heatmap.txt" matrix with image, (7.49 <= x && x < 12.55) ? floor(x+1.5)-0.5 : 0/0 ls 1

set cblabel "% AGU"
set cbtics 0.2 format "%.1f"
set xlabel "Distance from 3' end (codons)"
plot "data/ecoli_ser/AGT.heatmap.txt" matrix with image, (7.49 <= x && x < 12.55) ? floor(x+1.5)-0.5 : 0/0 ls 1

set cblabel "% AGC"
set cbtics 1 format "%.0f"
set xlabel " "
plot "data/ecoli_ser/AGC.heatmap.txt" matrix with image, (7.49 <= x && x < 12.55) ? floor(x+1.5)-0.5 : 0/0 ls 1


unset multiplot
