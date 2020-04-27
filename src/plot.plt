#! /usr/bin/env gnuplot

set encoding utf8
set tics font "Cantarell Bold, 50"
unset key
# set key box linewidth 0.1 width -4
# set grid
# set style line 1 \
#     linetype 3 linewidth .7 \
#     pointtype 2 pointsize 0.01

set output 'fig1.png'
set terminal png transparent \
    font 'Cantarell Bold' 50 \
    size 3840, 2160 \
    enhanced
set size 1.0, 0.8

set border linewidth 10

set xrange [0.0001:100]
set yrange [-50:70]

set ylabel "TEMP. ( DEG. C )"
set xlabel "TIME ( SECOND )"

set label "Curve A: a = 5 x 10^{-4} cm^{-1} deg^{-1}, b = 66.7 cm^{-1}, v_0 = 2000 cm sec^{-1} ;" at 0.001, 100
set label "Curve B: a = 2 x 10^{-4} cm^{-1} deg^{-1}, b = 667 cm^{-1}, v_0 = 200 cm sec^{-1} ;" at 0.001, 90
set label "Curve C: a = 5 x 10^{-4} cm^{-1} deg^{-1}, b = 6670 cm^{-1}, v_0 = 20 cm sec^{-1} . " at 0.001, 80

set xtics ("10^{-4}" 0.0001, "10^{-3}" 0.001, \
    "10^{-2}" 0.01, "10^{-1}" 0.1, "1" 1, "10^1" 10, "10^2" 100)
set logscale x

set label "A" at 0.001, 55
set label "B" at 0.001, 33
set label "C" at 0.001, 3
set label "T = T_0 + ( T_i - T_0 ) / [a b ( T_i - T_0 ) ln ( c t + 1 ) + 1]" \
    at 0.1, 55
set label "T_i = INITIAL TEMP." at 0.1, 48
set label "T_0 = AMBIENT TEMP." at 0.1, 41
plot 'fig1_1.txt' using 1:2 with lines linewidth 7 dt 1 lt 8, \
    'fig1_2.txt' using 1:2 with lines linewidth 7 dt 1 lt 8, \
    'fig1_3.txt' using 1:2 with lines linewidth 7 dt 1 lt 8

set output 'fig2.png'
set terminal png transparent \
    font 'Cantarell Bold' 14 \
    size 700, 1000 \
    enhanced
# set terminal pdf enhanced font 'Cantarell Bold' fontscale 0.2 size 1.4, 2.0
set xrange [230:340]
set yrange [100:200000]
set y2range [2.325:2.825]

set xtics 230 340 20
set mxtics 2
set ytics ("10^2" 100, "10^3" 1000, "10^4" 10000, "10^5" 100000)
set mytics 4

set xlabel "T ( K )"
set ylabel "e_s ( DYNE/CM^2 )"

set label "e_s = e_w ; T > 273 K" at 240, 150000 
set label "e_s = e_w + ( e_i - e_w ) ( 273 - T ) / 40 ; T < 273 K" \
    at 240, 120000

set label "e_i = SAT. VAP. PRESSURE OVER ICE" at 240, 90000 
set label "e_w = SAT. VAP. PRESSURE \n            OVER WATER" at 240, 75000

unset logscale x
set logscale y

plot 'fig2_1.txt' using 1:2 with lines linewidth 2.7 dt 1 lt 8

set output 'fig3.png'
set terminal png transparent \
    font 'Cantarell Bold' 14 \
    size 960, 560 \
    enhanced

set xrange [220:380]
set yrange [65:85]

set xtics 220 380 20
set mxtics 2
set ytics 65 85 5
set mytics 5

set y2tics ("2.4" 2.4, "2.5" 2.5, "2.6" 2.6, "2.7" 2.7, "2.8" 2.8)
set y2label "L(ERG/GM) x 10^{-10}"
set my2tics 4

set xlabel "T ( K )"
set ylabel "σ ( DYNE/CM )"

set label "σ" at 230, 81
set label "{/Sans:Bold L}" at 240, 83
set label "L = L_v ; T > 273 K" at 290, 84
set label "L = L_v + ( L_s - L_v ) ( 273 - T ) / 40 ; T < 273 K" at 290, 83
set label "L_v = HEAT OF VAPORIZATION" at 290, 81.9
set label "L_s = HEAT OF SUBLIMATION" at 290, 81
unset logscale y 
plot 'fig3_1.txt' using 1:2 with lines linewidth 2.5 dt 1 lt 8 axes x1y2, \
    'fig3_1.txt' using 1:3 with lines linewidth 2.5 dt 2 lt 8 axes x1y1

set output "fig4.png"
set terminal png transparent \
    font 'Cantarell Bold' 14 \
    size 430, 615 \
    enhanced

set xrange [230:340]
set yrange [50000:200000000]

set xtics 230 340 20
set mxtics 2
set ytics ("10^5" 100000, "10^6" 1000000, "10^7" 10000000, "10^8" 100000000)
set mytics 10

unset ylabel
unset y2label
unset y2tics
unset label

set label "u" at 240, 1700000
set label "v" at 260, 5750000
set label "u = L^2 M / K R T^2" at 290, 40000000
set label "v = R T / D M e_s" at 290, 30000000

set logscale y
plot 'fig4_1.txt' using 1:2 with lines linewidth 2.5 dt 1 lt 8, \
    'fig4_1.txt' using 1:3 with lines linewidth 2.5 dt 1 lt 8

set output "fig5.png"
set terminal png transparent \
    font 'Cantarell Bold' 50 \
    size 3840, 2160 \
    enhanced
    
set size 1.0, 0.85

set border linewidth 10
set tics font "Cantarell Bold, 50"

set xrange [0.0001:0.01]
set yrange [1:10]

set logscale x 10
set logscale y 10

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

plot 'fig5_1.txt' using 1:2 with lines lw 7 dt 5 lt 8, \
     'fig5_2.txt' using ($1*0.01):2 with lines lw 7 dt 1 lt 8, \
     'fig5_3.txt' using ($1*0.001):2 with lines lw 7 dt 1 lt 8
