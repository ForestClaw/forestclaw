
How to do the strong/weak scaling plots on Stampede:

1. Make sure that there is NO serialout=T in any configuration file.

2. Make a bunch of runs with nout=1 and an appropriate tfinal value using:
   #SBATCH -o swirl.%j.out
   #SBATCH -e swirl.%j.err
   I have a job script doc/stampede/swirl that's a good starting point.

3. Copy the .out files somewhere, call path/to/analyze.pl -s <LNUM>
   where LNUM is the number of cores for the smallest run.
   This creates six text files with the times.

4. Call path/to/plot.sh <LNUM> <HNUM> <MX> from the same diretory.
   This creates six eps files with the plots.
   If regrid_interval=0, the regrid plots will error and be empty.

For single-node MPI I have used LNUM=1, HNUM=16.
For multi-node MPI I have used LNUM=16, HNUM=256.
On subsequent runs, the core numbers should increase by powers of 2.

It's easy to tweak plot.sh to change gnuplot options, ranges, etc.

Important command line options for VTK output are --nout=something and
--prefix=$SCRATCH/some/path/$SLURM_JOB_NAME.$SLURM_JOB_ID --vtkout=3.
