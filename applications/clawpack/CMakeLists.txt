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
add_subdirectory(advection/2d/swirl)

## filament  (square, multi-block disk)
add_subdirectory(advection/2d/filament)

## annulus  (manifold, periodic in one direction)
add_subdirectory(advection/2d/annulus)

## latlong (manifold, periodic in one direction)
add_subdirectory(advection/2d/latlong)

## disk  (manifold, multi-block)
add_subdirectory(advection/2d/disk)

## torus  (manifold, periodic in both directions)
add_subdirectory(advection/2d/torus)

## replicated (square, multi-block)
add_subdirectory(advection/2d/replicated)

## hemisphere (manifold, multi-block)
add_subdirectory(advection/2d/hemisphere)

## sphere (cubed-sphere and pillow-disk)
add_subdirectory(advection/2d/sphere)

## ----------------------------------------------------------
## Other hyperbolic problems (acoustics, Euler, burgers, SWE)
## ----------------------------------------------------------

## Acoustics (on flat domains)
add_subdirectory(acoustics/2d/radial)
add_subdirectory(acoustics/2d/interface)

## Burgers
add_subdirectory(burgers/2d/pwconst)

## Shallow
add_subdirectory(shallow/2d/radialdam)
add_subdirectory(shallow/2d/bump)

## Euler
add_subdirectory(euler/2d/shockbubble)
add_subdirectory(euler/2d/quadrants)




## ---------------------Transport on the sphere ---------------------------- 
## Examples from suite described by Lauritzen, et al.  See:
##   "A standard test case suite for two-dimensional  linear
##    transport on the sphere: results from a collection 
##    of state-of-the-art schemes", Lauritzen, et al.
##    Geosciences Model Development, 2014. 
##    http://www.geosci-model-dev.net/7/105/2014/gmd-7-105-2014.html
## -------------------------------------------------------------------------

add_subdirectory(transport/gaussian)
#add_subdirectory(transport/slotted_disk)
#add_subdirectory(transport/correlatedcb)


## ---------------------Transport on the sphere ---------------------------- 
## Example : Don't using the 'app' for configuring options;  
##           This is how you might call ForestClaw from another
##           program.
## -------------------------------------------------------------------------
#add_subdirectory(advection/2d/torthem)


