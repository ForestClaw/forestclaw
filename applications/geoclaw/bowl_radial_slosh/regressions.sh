#!/bin/sh
# absolute path to application we are testing
application=$FCLAW_APPLICATIONS_BUILD_DIR/geoclaw/bowl_radial_slosh/bowl_radial_slosh

# change to source dir for working directory
cd $FCLAW_APPLICATIONS_SRC_DIR/geoclaw/bowl_radial_slosh

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini  --radial:regression-check=radial_regression_values_0.csv --slosh:regression-check=slosh_regression_values_0.csv|| exit 1
