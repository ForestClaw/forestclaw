# --------------------------------------------
# ThunderEgg solver
# --------------------------------------------

## Elliptic example
add_subdirectory(elliptic/2d/poisson)
add_subdirectory(elliptic/3d/poisson)

## Heat equation
add_subdirectory(elliptic/2d/heat)

## Allencahn equation
add_subdirectory(elliptic/2d/allencahn)

## Crystal growth
add_subdirectory(elliptic/2d/phasefield)

## Heat Phasefield two solver example
add_subdirectory(elliptic/2d/heat_phasefield)