#!/bin/sh
# absolute path to application we are testing
application=$(pwd)/applications/clawpack/advection/2d/hemisphere/hemisphere

# change to source dir for working directory
cd $srcdir/applications/clawpack/advection/2d/hemisphere/

# run programs, exit script with nonzero on failure
mpirun -n 2 $application -F regression.ini --user:claw-version=4 --user:example=0 || exit 1
mpirun -n 2 $application -F regression.ini --user:claw-version=5 --user:example=0 || exit 1
mpirun -n 2 $application -F regression.ini --user:claw-version=4 --user:example=1 || exit 1
mpirun -n 2 $application -F regression.ini --user:claw-version=5 --user:example=1 || exit 1