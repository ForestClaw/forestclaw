# ForestClaw

http://www.forestclaw.org


ForestClaw (http://www.forestclaw.org) is a parallel, multi-block adaptive finite volume
code  for solving PDEs on a hierarchy of logically Cartesian meshes.  

Features of ForestClaw : 

* A multi-resolution hierarchy of grids is stored  as a composite structure of non-
overlapping fixed sized  grids which are stored as leaves in a forest of quad- or octrees.

* Based on the highly scalable meshing library p4est (http://www.p4est.org)
 
* Solvers for hyperbolic PDEs are available, including Clawpack 4.x, Clawpack 5.0 and
GeoClaw (see http://www.clawpack.org).
    
* Users can easily extend ForestClaw with their own solvers. 
    
* Fully parallel, using MPI distributed memory parallelism
    
* Visualization tools include Matlab scripts, VisClaw (from Clawpack) and VTK output.

For installation instructions, please visit our GitHub Wiki at https://github.com/ForestClaw/ForestClaw/wiki.

For more information on ForestClaw, visit our website at http://www.forestclaw.org.

ForestClaw is free software under a BSD-style license (see COPYING). The p4est library is
released under the GPL 2.0 (or later).  As it is generally linked with the ForestClaw
code,  binary distribution falls under GPL as well.

---
