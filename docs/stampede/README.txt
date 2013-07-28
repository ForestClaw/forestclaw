
How to do the strong/weak scaling plots:

1. Make sure that there is NO serialout=T in any configuration file.

2. Make a bunch of runs with nout=1 and an appropriate tfinal value.

3. Copy the .out files somewhere, call path/to/analyze.pl -s <LNUM>
   where CNUM is the number of cores for the smallest run.
   This creates six text files with the times.

4. Call path/to/plot.sh <LNUM> <HNUM> <MX> from the same diretory.
   This creates six eps files with the plots.

For single-node MPI I have used LNUM=1, HNUM=16.
For multi-node MPI I have used LNUM=16, HNUM=256.

It's easy to tweak plot.sh to change gnuplot options, ranges, etc.
