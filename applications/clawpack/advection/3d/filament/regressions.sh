#!/bin/sh
# absolute path to application we are testing
application=$FCLAW_APPLICATIONS_BUILD_DIR/clawpack/advection/3d/filament/filament

# change to source dir for working directory
cd $FCLAW_APPLICATIONS_SRC_DIR/clawpack/advection/3d/filament/

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression.ini || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:example=0 || exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:example=1|| exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:example=2|| exit 1
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application -F regression_map.ini --user:example=3|| exit 1

