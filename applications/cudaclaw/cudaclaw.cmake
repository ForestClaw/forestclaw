# ----------------------------------
# CudaClaw library and examples
# ----------------------------------


# Scalar advection
add_subdirectory(cudaclaw/advection/2d/swirl)

# Shallow water wave equations
add_subdirectory(cudaclaw/shallow/2d/bump)
#add_subdirectory(cudaclaw/shallow/2d/radialdam)
#add_subdirectory(cudaclaw/shallow/2d/tsunami)

# Acoustics
add_subdirectory(cudaclaw/acoustics/2d/radial)

# Euler equations
add_subdirectory(cudaclaw/euler/2d/shockbubble)
