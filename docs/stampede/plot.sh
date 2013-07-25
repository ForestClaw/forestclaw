#! /bin/bash

if test -z "$3" ; then
	echo "Arguments: startnumprocs, maxnumprocs, mxy"
	exit 1
fi
SNP="$1"
MNP="$2"
MXY="$3"

(cat <<EOF
# This is the Gnuplot script

set term postscript color solid
set style data linespoints

set logscale xy
set xlabel "Number of cores"
set format y "%.2e"

# Strong scaling plots

set xtics scale 0,0 $SNP,2,$MNP
set title "Advection, strong scaling for mx = my = $MXY, uniform grid"

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

# Weak scaling plots

set xtics scale 0,0 $SNP,4,$MNP
set title "Advection, weak scaling for mx = my = $MXY, uniform grid"

set ylabel "Compute time per time step (seconds)"
set output "advance_weak.eps"
plot "advance_weak.txt" using 1:2 title "Levels 1-3", \
     "advance_weak.txt" using 1:3 title "Levels 2-4", \
     "advance_weak.txt" using 1:4 title "Levels 3-5", \
     "advance_weak.txt" using 1:5 title "Levels 4-6", \
     "advance_weak.txt" using 1:6 title "Levels 5-7", \
     "advance_weak.txt" using 1:7 title "Levels 6-8"

set ylabel "Exchange time per time step (seconds)"
set output "exchange_weak.eps"
plot "exchange_weak.txt" using 1:2 title "Levels 1-3", \
     "exchange_weak.txt" using 1:3 title "Levels 2-4", \
     "exchange_weak.txt" using 1:4 title "Levels 3-5", \
     "exchange_weak.txt" using 1:5 title "Levels 4-6", \
     "exchange_weak.txt" using 1:6 title "Levels 5-7", \
     "exchange_weak.txt" using 1:7 title "Levels 6-8"

EOF
) | gnuplot
