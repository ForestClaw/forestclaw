#! /bin/sh

MINL=1
MAXL=5

NODES="1 2 4"


for NODE in $NODES ; do
	# This defines the processor count
	echo "# nodes $NODE"
        perl -e "print 16 * $NODE"

for LEVEL in `seq $MINL $MAXL` ; do
	# This defines the problem size
	#echo "# level $LEVEL"

	# not computing the level facter
        grep "Max/P" "mx32_omp_16_level_${LEVEL}_mpi_${NODE}.out" | \
		cut -d ' ' -f 5,6 | awk "{ printf (\" %g\", \$2 / \$1 / (1**$LEVEL)); }"
done
	echo
done
