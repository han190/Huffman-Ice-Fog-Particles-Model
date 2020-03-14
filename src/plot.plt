#! /usr/bin/env gnuplot

set terminal pdf enhanced font 'Cantarell Bold' \
    fontscale 0.236029 size 2.00in, 1.236in
set key box linewidth 0.1
set grid
set style line 1 \
    linetype 1 linewidth 1 \
    pointtype 7 pointsize 0.01

set output 'fig.pdf'
set logscale x
plot 'fig1_1.txt' using 1:2 with linespoints linestyle 1 lt rgb 'red', \
     'fig1_2.txt' using 1:2 with linespoints linestyle 1 lt rgb 'blue', \
     'fig1_3.txt' using 1:2 with linespoints linestyle 1 lt rgb 'green'

filename = 'fig2_1.txt'
# set logscale xy
plot filename using 3:5 with linespoints linestyle 1

plot filename using 4:7 with linespoints linestyle 1

