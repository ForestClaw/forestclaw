## ------------------------- Clawpack examples ----------------------
##
## Hyperbolic problems using Clawpack 4.6 and 5.0 
## 
## Many of these examples have been taken from the original 
## Clawpack software package, developed by R. J. LeVeque, 
## and many others. (www.clawpack.org)
## 
## For descriptions of the wave propagation algorithm used 
## in Clawpack, see : 
## 
## See : "Finite volume methods for hyperbolic problems", 
##        R. J. LeVeque, Cambridge University Press, 2002.
## -----------------------------------------------------------------

## ----------------
## Scalar advection 
## ----------------

## swirl  (single block, square)
add_subdirectory(clawpack/advection/2d/swirl)

## filament  (square, multi-block disk)
add_subdirectory(clawpack/advection/2d/filament)

## filament_swirl
add_subdirectory(clawpack/advection/2d/filament_swirl)

## annulus  (manifold, periodic in one direction)
add_subdirectory(clawpack/advection/2d/annulus)

## latlong (manifold, periodic in one direction)
add_subdirectory(clawpack/advection/2d/latlong)

## disk  (manifold, multi-block)
add_subdirectory(clawpack/advection/2d/disk)

## torus  (manifold, periodic in both directions)
add_subdirectory(clawpack/advection/2d/torus)

## replicated (square, multi-block)
add_subdirectory(clawpack/advection/2d/replicated)

## hemisphere (manifold, multi-block)
add_subdirectory(clawpack/advection/2d/hemisphere)

## sphere (cubed-sphere and pillow-disk)
add_subdirectory(clawpack/advection/2d/sphere)

## periodic (square,periodic in both directions, constant velocity)
add_subdirectory(clawpack/advection/2d/periodic)

## --------------------- Advection on a sphere ---------------------------- 
## Examples from suite described by Lauritzen, et al.  See:
##   "A standard test case suite for two-dimensional  linear
##    transport on the sphere: results from a collection 
##    of state-of-the-art schemes", Lauritzen, et al.
##    Geosciences Model Development, 2014. 
##    http://www.geosci-model-dev.net/7/105/2014/gmd-7-105-2014.html
##
## Like the examples above, these three examples all solve the color
## color equation
## -------------------------------------------------------------------------

add_subdirectory(clawpack/advection/2d/gaussian)
add_subdirectory(clawpack/advection/2d/correlatedcb)
add_subdirectory(clawpack/advection/2d/slotted_disk)

## ----------------------------- Transport --------------------------------- 
## Example : Solve the variable velocity transport equation
##
##                   q_t + div(u q) = 0
##
## These examples use the fwave approach and can be used to test
## conservation.
##
## -------------------------------------------------------------------------

add_subdirectory(clawpack/transport/2d/sphere)
add_subdirectory(clawpack/transport/2d/torus)
add_subdirectory(clawpack/transport/2d/square)

## -------------------------------- Rays --------------------------------- 
## Test example including rays
## -------------------------------------------------------------------------

add_subdirectory(clawpack/advection/2d/swirl_rays)

## ----------------------------------------------------------
## Other hyperbolic problems (acoustics, Euler, burgers, SWE)
## ----------------------------------------------------------

## Acoustics (on flat domains)
add_subdirectory(clawpack/acoustics/2d/radial)
add_subdirectory(clawpack/acoustics/2d/interface)

## Burgers
add_subdirectory(clawpack/burgers/2d/pwconst)

## Shallow
add_subdirectory(clawpack/shallow/2d/radialdam)
add_subdirectory(clawpack/shallow/2d/bump)

## Euler
add_subdirectory(clawpack/euler/2d/shockbubble)
add_subdirectory(clawpack/euler/2d/quadrants)
add_subdirectory(clawpack/euler/2d/triple)

## --------------------- Miscellaneous ---------------------------- 
## Example : Don't using the 'app' for configuring options;  
##           This is how you might call ForestClaw from another
##           program.
## -------------------------------------------------------------------------
#add_subdirectory(clawpack/advection/2d/torthem)


