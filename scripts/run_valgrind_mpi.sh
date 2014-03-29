#!/bin/tcsh

# run_valgrind.sh 2 swirl -F fclaw2d_user.ini
set nprocs = $1
set target = $2
set flags = "$3 $4 "

echo "${target} ${flags}"


set supp_c = --suppressions=${FORESTCLAW}/fclaw2d_vg.supp
set gen_supp = --gen-suppressions=all
set log_file = --log-file=${target}.out
set track = --track-origins=yes
set leak = --leak-check=full
set reach = --show-reachable=yes
set debug = "--vgdb-error=0 --db-attach=yes"

set supp_user = --suppressions=user.supp

set ops = "${leak} ${reach} ${track} ${supp_c} ${gen_supp}"

dsymutil ${target}
valgrind -v  ${ops} mpirun -np $nprocs ${target} ${flags}
