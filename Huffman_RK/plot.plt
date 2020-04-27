#! /usr/bin/env gnuplot

set encoding utf8
unset key

set output "fig.png"
set terminal png transparent \
    font 'Cantarell Bold' 50 \
    size 3840, 2160 \
    enhanced
    
set size 1.0, 0.85

set border linewidth 10
set tics font "Cantarell Bold, 50"

set xrange [0.0001:0.01]
set yrange [1:10]

set xlabel "TIME ( SECONDS )" offset 0, 1
set ylabel "SATURATION RATIO"

set xtics mirror
set xtics (" " 0.0001, " " 0.0002, " " 0.0003, " " 0.0004, \
    " " 0.0005, " " 0.0006, " " 0.0007, " " 0.0008, " " 0.0009, \
    " " 0.001, " " 0.002, " " 0.003, " " 0.004, " " 0.005, \
    " " 0.006, " " 0.007, " " 0.008, " " 0.009, " " 0.01)
# set mxtics 10
set ytics ("1" 1, "2" 2, "3" 3, "4" 4, "5" 5, "6" 6, \
    "7" 7, "8" 8, "9" 9, "10" 10)
    
set label "A" at 0.0013, 2.8
set label "B" at 0.001, 1.4
set label "C" at 0.002, 4.25

set label "A:10^{-4}" at 0.0001, 11 left
set label "10^{-3}" at 0.001, 11 left
set label "10^{-2}" at 0.01, 11 left
set label "B:10^{-2}" at 0.0001, 12.5 left
set label "10^{-1}" at 0.001, 12.5 left
set label "10^{0}" at 0.01, 12.5 left
set label "C:10^{-1}" at 0.0001, 14.2 left
set label "10^{0}" at 0.001, 14.2 left
set label "10^{1}" at 0.01, 14.2 left
    
set logscale x 10
set logscale y 10
    
# plot \
#     'fig1.txt' using 1:2 with lines linewidth 7 dt 1 lt 8, \
#     'fig2.txt' using ($1*0.01):2 with lines linewidth 7 dt 1 lt 8, \
#     'fig3.txt' using ($1*0.001):2 with lines linewidth 7 dt 1 lt 8

plot \
    'fig1.txt' using 1:2 with points pointtype 7 pointsize .05, \
    'fig2.txt' using ($1*0.01):2 with points pointtype 7 pointsize .05, \
    'fig3.txt' using ($1*0.001):2 with points pointtype 7 pointsize .05
