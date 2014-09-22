#!/bin/tcsh

set target = $1
set flags = "$2 $3"

echo "${target} ${flags}"


set supp_c = --suppressions=$1.supp
set gen_supp = --gen-suppressions=all
set log_file = --log-file=${target}.vglog
set track = --track-origins=yes
set leak = --leak-check=full
set reach = --show-reachable=yes
set debug = "--vgdb-error=0 --db-attach=yes"

set supp_user = --suppressions=user.supp

set ops = "${leak} ${reach} ${track} ${supp_c} ${gen_supp}"

dsymutil ${target}
valgrind -v  ${ops} ${target} ${flags}
