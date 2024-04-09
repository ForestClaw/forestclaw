
1.  Get a unit problem (exponential spike problem, for example) in a [0,1]x[0,1]
domain.

2.  Refine the right hand between levels 3 and 5 (for example) so that on a single block
we have roughly 4000 patches, and the unit problem can be solved in about a second (or less)

3.  Run the unit problem on 1, 4, 16, 64, 256 procs (or possibly 1,2,4,8,16,32) 

4.  Then,duplicate the unit problem on a 2x2 grid of blocks and run on 4,16,32,64,256 procs

5.  Continue for block grids 4x4, 8x8 and 16,16

6.  Produces tables like the one below with various metrics : 


Table 1 : Grids per processor (show below)

Table 2 : Wall clock time

Table 3 : Cells/time 

Table 4 : Efficiency  actual speed-up/expected speed-up.


Plots : Strong scaling (down the columns).  Weak scaling is along the diagonals (constant work per processor)




—————————————————————————————————————————————————————
 Procs |      1x1      2x2      4x4      8x8    16x16
—————————————————————————————————————————————————————
    1  |     4096      ---      ---      ---      ---
    4  |     1024     4096      ---      ---      ---
   16  |      256     1024     4096      ---      ---
   64  |       64      256     1024     4096      ---
  256  |       16       64      256     1024     4096
—————————————————————————————————————————————————————


BC conditions
-------------


Dirichlet : 

u = 7 on teh boundary :

(Ui + Ug)/2 = 7   Ui = u interior;  Ug = u ghost value

Ug = 2*7 - Ui


The user should the the ghost value to Ug.  They do this for each ghost cell on the patch that 
is passed into the BC routine.  The interior of this patch is set to zero.

Then you apply the Laplacian to the interior cells in the patch.  Store the result in 
the interior cells at the boundary.    Add values for corner cells.

The result patch of values should be subtracted from the right hand side supplied by the user to 
construct a new right hand side.

Solve the resulting and Voila, the solution should have the right BCs!




