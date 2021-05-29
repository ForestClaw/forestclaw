# ----------------------------------
# CudaClaw library and examples
# ----------------------------------


# Scalar advection
add_subdirectory(clawpack/advection/2d/swirl_cuda)

# Shallow water wave equations
#include applications/clawpack/shallow/2d/bump_cuda/Makefile.am
##include applications/clawpack/shallow/2d/radialdam_cuda/Makefile.am
##include applications/clawpack/shallow/2d/tsunami_cuda/Makefile.am

# Acoustics
##include applications/clawpack/acoustics/2d/radial_cuda/Makefile.am

# Euler equations
##include applications/clawpack/euler/2d/shockbubble_cuda/Makefile.am
