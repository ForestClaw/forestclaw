------------------
Phasefield Solver
------------------

This code solves the isotropic phasefield equations (1) and (2) in Murray et al. (mu-wh-gl.pdf).  

Two main routines are supplied : phase_cc.f and phase_iterative.f90. 

Edit TARGET in the Makefile to select version.

-----------------------
Version 1:  phase_cc.f 
-----------------------
Solves equation (1) followed by (2) for u^{n+1} and phi^{n+1}.   In (1), dphi/dt is approximated as : 

          dphi/dt = (phi^{n} - phi^{n-1})/dt

After solving (1), u^{n+1} is then used in right hand of (2). We are required to store the previous time level phi^{n-1}.

Each time step requires two solves (one for each equation (1) and (2)).  Each solve is a call to FISHPACK hstcrt.   

Parameters for this version can be set in `phase.dat'. 

------------------------------
Version 2:  phase_iterative.f  
------------------------------
Solves a block 2x2 system (1) and (2) iteratively.  The main difference between Version 1 and Version 2 is how dphi/dt is computed.  In version 2, we compute dphi/dt implicitly : 

       dphi/dt = (phi^{n+1} - phi^{n})/dt

This form is what is needed for an AMR mesh, since we won't typically have access to previous time levels.  A second, minor difference is that the value of u on the right hand side of the phi equation is treated implicitly (for Gauss-Seidel) or explicitly (Jacobi).  

This version converges with both Jacobi or Gauss-Seidel.  Both however, are quite slow and require about 50 (Jacobi) or 20 (Gauss-Seidel) iterations per time step. Each iteration requires two calls to hstcrt.f.  

See phasefield_equations.{tex,pdf} for the block form of the equations.  


Parameters for this version can be set in `phase_iterative.dat'. 


-------------------------------
FISHPACK
-------------------------------

This phasefield solver uses the cell-centered Poisson solver hstcrt.f, available in the Fishpack library.  Fishpack can be downloaded here : 

    https://www2.cisl.ucar.edu/resources/legacy/fishpack

Fishpack is an ancient code that defaults to single precision. To restore full double precision, compile Fishpack with gfortran flags : 

Linux:
------
    F90 := gfortran -fdefault-real-8 -O2
    CPP := gfortran -cpp

This can be set, along with other customizable settings, in 'make.inc'.  


-------------------------------
Visualization
-------------------------------

Output is written to `fort.qXXXX` files and `fort.tXXXX'.  The output can be visualized in Matlab using VisClaw (download as part of the clawpack package www.clawpack.org).  Also required are Matlab routines in forestclaw.  (www.forestclaw.org).  

***********************************************************************






