#!/bin/sh
# absolute path to application we are testing
application=$(pwd)/applications/clawpack/advection/2d/hemisphere/hemisphere

# change to source dir for working directory
cd $srcdir/applications/clawpack/advection/2d/hemisphere/

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini --user:claw-version=4 --user:example=0 || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini --user:claw-version=5 --user:example=0 || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini --user:claw-version=4 --user:example=1 || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini --user:claw-version=5 --user:example=1 || exit 1


echo $FCLAW_MPI_RUN $FCLAW_MPI_TEST_FLAGS
exit 1