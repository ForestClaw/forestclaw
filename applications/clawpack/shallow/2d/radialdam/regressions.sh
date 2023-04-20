#!/bin/sh
# absolute path to application we are testing
application=$FCLAW_APPLICATIONS_BUILD_DIR/clawpack/shallow/2d/radialdam/radialdam

# change to source dir for working directory
cd $FCLAW_APPLICATIONS_SRC_DIR/clawpack/shallow/2d/radialdam/

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini     --user:claw-version=4 --user:example=0 --trapfpe=F  --regression-check=regression_values_0.csv || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini     --user:claw-version=5 --user:example=0 --trapfpe=F  --regression-check=regression_values_1.csv || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:claw-version=4 --user:example=1 --trapfpe=F  --regression-check=regression_values_2.csv || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:claw-version=5 --user:example=1 --trapfpe=F  --regression-check=regression_values_3.csv || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:claw-version=4 --user:example=2 --trapfpe=F  --regression-check=regression_values_4.csv || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:claw-version=5 --user:example=2 --trapfpe=F  --regression-check=regression_values_5.csv || exit 1
