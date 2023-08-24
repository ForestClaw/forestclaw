#!/bin/sh
# absolute path to application we are testing
application=$FCLAW_APPLICATIONS_BUILD_DIR/clawpack/advection/2d/annulus/annulus

# change to source dir for working directory
cd $FCLAW_APPLICATIONS_SRC_DIR/clawpack/advection/2d/annulus/

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini --user:claw-version=4  --regression-check=regression_values_0.csv || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini --user:claw-version=5  --regression-check=regression_values_1.csv || exit 1
