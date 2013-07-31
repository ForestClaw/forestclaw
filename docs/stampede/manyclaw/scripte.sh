#! /bin/sh

MINL=1
MAXL=8

NODES="1 2 4 8"


for NODE in $NODES ; do
	# This defines the processor count
	echo "# nodes $NODE"
        perl -e "print 16 * $NODE"

for LEVEL in `seq $MINL $MAXL` ; do
	# This defines the problem size
	#echo "# level $LEVEL"

	# Not computing the level factor
        grep "Max/P" "mx32_omp_16_level_${LEVEL}_mpi_${NODE}.out" | \
		cut -d ' ' -f 8,9 | awk "{ printf (\" %g\", \$2 / \$1 / (1**$LEVEL)); }"
done
	echo
done
