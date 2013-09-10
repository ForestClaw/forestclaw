#! /bin/sh

cat *.out | grep Procs | sort -t ' ' -k 3 -n | \
awk '/Procs/ { \
	t=$6/$5; if (tf == 0) { tf = t * $3; } \
	printf ("%8d %8.3g %8.1f %8.1f\n", $3, t, t*$3, t*$3 / tf * 100.); }' \
| tee strong.txt
