# ----------------------------------
# CudaClaw library and examples
# ----------------------------------


# Scalar advection
add_subdirectory(clawpack/advection/2d/swirl_cuda)

# Shallow water wave equations
add_subdirectory(clawpack/shallow/2d/bump_cuda)
add_subdirectory(clawpack/shallow/2d/radialdam_cuda)
add_subdirectory(clawpack/shallow/2d/tsunami_cuda)

# Acoustics
add_subdirectory(clawpack/acoustics/2d/radial_cuda)

# Euler equations
add_subdirectory(clawpack/euler/2d/shockbubble_cuda)
