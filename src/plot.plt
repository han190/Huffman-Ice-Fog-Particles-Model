#! /usr/bin/env gnuplot

set encoding utf8
unset key
set output "saturation.png"
set terminal png transparent \
    font 'Cantarell Bold' 50 \
    size 3840, 2160 \
    enhanced
set size 1.0, 0.85
set border linewidth 10
set tics font "Cantarell Bold, 50"
set xrange [0.0001:0.01]
set yrange [1:10]
set xlabel "TIME (SECOND)" offset 0, 1
set ylabel "SATURATION RATIO"
set xtics mirror
set xtics (" " 0.0001, " " 0.0002, " " 0.0003, " " 0.0004, \
    " " 0.0005, " " 0.0006, " " 0.0007, " " 0.0008, " " 0.0009, \
    " " 0.001, " " 0.002, " " 0.003, " " 0.004, " " 0.005, \
    " " 0.006, " " 0.007, " " 0.008, " " 0.009, " " 0.01)
set ytics ("1" 1, "2" 2, "3" 3, "4" 4, "5" 5, "6" 6, \
    "7" 7, "8" 8, "9" 9, "10" 10)
    
set label "C" at 0.003, 7.0
set label "B" at 0.002, 3.1
set label "A" at 0.003, 2.8
set label "A:10^{-4}" at 0.0001, 11 left
set label "10^{-3}" at 0.001, 11 center
set label "10^{-2}" at 0.01, 11 center
set label "B:10^{-2}" at 0.0001, 12.5 left
set label "10^{-1}" at 0.001, 12.5 center
set label "10^{0}" at 0.01, 12.5 center
set label "C:10^{-1}" at 0.0001, 14.2 left
set label "10^{0}" at 0.001, 14.2 center
set label "10^{1}" at 0.01, 14.2 center
set logscale x 10
set logscale y 10
    
plot \
    'saturation_A.dat' using 1:2 with lines \
    linewidth 7 dt 1 lt 8 title "A", \
    'saturation_B.dat' using ($1*0.01):2 with lines \
    linewidth 7 dt 1 lt 8 title "B", \
    'saturation_C.dat' using ($1*0.001):2 with lines \
    linewidth 7 dt 1 lt 8 title "C"

reset

set encoding utf8
unset key
set output 'temperature.png'
set terminal png transparent \
    font 'Cantarell Bold' 50 \
    size 3840, 2160 \
    enhanced
set size 1.0, 0.8
set border linewidth 10
set xrange [0.0001:100]
set yrange [-50:70]

set ylabel "TEMP. (DEG. C)"
set xlabel "TIME (SECOND)"
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

plot \
    'temp_A.dat' using 1:2 with lines linewidth 7 dt 1 lt 8, \
    'temp_B.dat' using 1:2 with lines linewidth 7 dt 1 lt 8, \
    'temp_C.dat' using 1:2 with lines linewidth 7 dt 1 lt 8


