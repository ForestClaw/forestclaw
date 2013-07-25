#! /bin/bash

(cat <<EOF
# This is the Gnuplot script

set term postscript color solid
set style data linespoints

set logscale xy
set xlabel "Number of cores"
set xtics scale 0,0 1,2,16
set format y "%.2e"

set ylabel "Compute time per time step (seconds)"
set output "advance_strong.eps"
plot "advance_strong.txt" using 1:2 title "Level 1", \
     "advance_strong.txt" using 1:3 title "Level 2", \
     "advance_strong.txt" using 1:4 title "Level 3", \
     "advance_strong.txt" using 1:5 title "Level 4", \
     "advance_strong.txt" using 1:6 title "Level 5", \
     "advance_strong.txt" using 1:7 title "Level 6", \
     "advance_strong.txt" using 1:8 title "Level 7", \
     "advance_strong.txt" using 1:9 title "Level 8"

set ylabel "Exchange time per time step (seconds)"
set output "exchange_strong.eps"
plot "exchange_strong.txt" using 1:2 title "Level 1", \
     "exchange_strong.txt" using 1:3 title "Level 2", \
     "exchange_strong.txt" using 1:4 title "Level 3", \
     "exchange_strong.txt" using 1:5 title "Level 4", \
     "exchange_strong.txt" using 1:6 title "Level 5", \
     "exchange_strong.txt" using 1:7 title "Level 6", \
     "exchange_strong.txt" using 1:8 title "Level 7", \
     "exchange_strong.txt" using 1:9 title "Level 8"

EOF
) | gnuplot
