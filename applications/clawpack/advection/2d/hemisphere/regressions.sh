#!/bin/sh
# absolute path to application we are testing
application=$(pwd)/applications/clawpack/advection/2d/hemisphere/hemisphere

# change to source dir for working directory
cd $srcdir/applications/clawpack/advection/2d/hemisphere/

# get runner
if command -v mpirun
then
	runner="mpirun -n 2"
fi

# run programs, exit script with nonzero on failure
$runner $application -F regression.ini --user:claw-version=4 --user:example=0 || exit 1
$runner $application -F regression.ini --user:claw-version=5 --user:example=0 || exit 1
$runner $application -F regression.ini --user:claw-version=4 --user:example=1 || exit 1
$runner $application -F regression.ini --user:claw-version=5 --user:example=1 || exit 1