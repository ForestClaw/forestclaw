#!/bin/bash

# Exit on any command failure
set -e

# doesn't have indirect cases, but checks that all calls are made correctly
mpirun -n 2 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on square with refinement"
mpirun -n 2 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on square with refinement coarse interior"
mpirun -n 2 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement on one block"
mpirun -n 2 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement on all but one block"
mpirun -n 2 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement coarse interior"
mpirun -n 2 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement 2"

# one patch per proc, each of these will have indirect cases
mpirun --oversubscribe -n 31 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on square with refinement"
mpirun --oversubscribe -n 52 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on square with refinement coarse interior"
mpirun --oversubscribe -n 28 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement on one block"
mpirun --oversubscribe -n 52 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement on all but one block"
mpirun --oversubscribe -n 61 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement coarse interior"
mpirun --oversubscribe -n 25 src/patches/clawpatch/clawpatch.TEST --indirect --test-case="2d clawpatch ghost fill on 2x2 brick with refinement 2"

