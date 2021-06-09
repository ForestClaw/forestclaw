## ----------------------- Under construction -----------------------------
##  These are not guaranteed to compile or may not be available with 
## standard distribution
## ------------------------------------------------------------------------

## Shallow water on the sphere (D. Calhoun)
#add_subdirectory(clawpack/shallow/2d/sphere)

## Simple examples to test options handling, etc, but no solvers.
#add_subdirectory(dummy)

## transport (square)
#add_subdirectory(clawpack/transport/2d/swirl)

## swirl_basilisk (square)
#add_subdirectory(clawpack/advection/2d/swirl_basilisk)

## torus_error (original)
#add_subdirectory(clawpack/advection/2d/torus_error)

## torus (original)
#add_subdirectory(clawpack/advection/2d/torus)

## periodic (original)
#add_subdirectory(clawpack/advection/2d/periodic)

## bump (original)
add_subdirectory(clawpack/shallow/2d/bump)

## Conservative advection (transport) on a square
#add_subdirectory(clawpack/transport/2d/square)

## tsunami (SGN)
#add_subdirectory(clawpack/shallow/2d/tsunami)

## Shallow (tsunami, used for SGN equations)
#add_subdirectory(clawpack/shallow/2d/tsunami)

## -----------------------------------------------------------
## Some third order examples
## -----------------------------------------------------------

## Burgers
#add_subdirectory(clawpack/burgers/2d/burgers_order3)

## Advection 
#add_subdirectory(clawpack/advection/2d/adv_order3)


## ---------------------------
## No solver
## ---------------------------
## This needs work because all exchange functions (copy, average, 
## interpolate) are stored in solvers, and not in a patch.  To 
## fix this, the patch needs to know something about data layout.  
## The clawpatch currently doesn't know about data layout. 

# add_subdirectory(no_solver)

## -------------------------------------------------
## Checks for metric terms (curvature, normals, etc)
## --------------------------------------------------
#add_subdirectory(metric/2d/all_mappings)
#add_subdirectory(metric/2d/mesh)



