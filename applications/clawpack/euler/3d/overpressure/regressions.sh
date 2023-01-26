#!/bin/sh
# absolute path to application we are testing
application=$FCLAW_APPLICATIONS_BUILD_DIR/clawpack/euler/3d/overpressure/overpressure

# change to source dir for working directory
cd $FCLAW_APPLICATIONS_SRC_DIR/clawpack/euler/3d/overpressure/

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)

$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application mpirun -n 4 overpressure -F cart.ini --tfinal=0.02 --nout=1 --output=F

$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application mpirun -n 4 overpressure -F cart_bump.ini --tfinal=0.02 --nout=1 --output=F

$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application mpirun -n 4 overpressure -F nomap.ini --tfinal=0.02 --nout=1 --output=F

$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $application mpirun -n 4 overpressure -F latlong.ini --tfinal=0.01 --nout=1 --output=F --clawpatch:mz=32

