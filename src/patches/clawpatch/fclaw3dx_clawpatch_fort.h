/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FCLAW3DX_CLAWPATCH_FORT_H
#define FCLAW3DX_CLAWPATCH_FORT_H

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw_global;
struct fclaw_patch;

struct fclaw_patch_transform_data;  /* Should be replaced by long int?  */

/**
 * @file 
 * Typedefs for clawpatch fortran functions
 * 
 * Functions defined here are implemented in individual solvers (clawpack 4.6 and 
 * clawpack 5.0) 
 */

#if 0
/* Fix syntax highlighting */
#endif

/** @{ @name Ghost filling - patch specific *//*----------------------------------------*/

/**
 * @brief Copies ghost data from a face-neighboring grid on the same level
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in,out] qthis the solution of this patch
 * @param[in]     qneighbor the solution of the neighbor patch
 * @param[in]     iface the interface that the neighbor is on
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_copy_face_t)(const int* mx, 
                                                   const int* my, 
                                                   const int* mz,
                                                   const int* mbc, 
                                                   const int* meqn,
                                                   double qthis[],
                                                   double qneighbor[], 
                                                   const int* iface,
                                                   struct fclaw_patch_transform_data** transform_ptr);

/**
 * @brief Averages values from a face-neighboring fine grid
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in,out] qcoarse the solution of this patch
 * @param[in]     qfine the solution of the fine neighbor patch
 * @param[in]     areacoarse the area of cells in this patch
 * @param[in]     areafine the area of cells in the fine neighbor patch
 * @param[in]     idir Face orientation - 0 for x-faces; 1 for y-faces [0-1]
 * @param[in]     iside the interface of the fine neighbor patch
 * @param[in]     num_neighbors the number of neighbors
 * @param[in]     refratio the refinement ratio
 * @param[in]     igrid the index of the fine neighbor in the child array
 * @param[in]     manifold true if using mainifold
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_average_face_t)(const int* mx, 
                                                      const int* my, 
                                                      const int* mz,
                                                      const int* mbc,    
                                                      const int* meqn,
                                                      double qcoarse[],
                                                      double qfine[],
                                                      double areacoarse[], 
                                                      double areafine[],
                                                      const int* idir, 
                                                      const int* iside,
                                                      const int* num_neighbors,
                                                      const int* refratio, 
                                                      const int* igrid,
                                                      const int* manifold, 
                                                      struct fclaw_patch_transform_data** transform_ptr);
/**
 * @brief Interpolates values from a face-neighboring coarse grid
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in]     qcoarse the solution of the coarse neighbor patch
 * @param[in,out] qfine the solution of this patch
 * @param[in]     areacoarse the area of cells in the coarse neighbor patch
 * @param[in]     arefine the area of cells in the this patch
 * @param[in]     idir Face orientation - 0 for x-faces; 1 for y-faces [0-1]
 * @param[in]     iface_coarse the interface of the coarse neighbor patch
 * @param[in]     num_neighbors the number of neighbors
 * @param[in]     refratio the refinement ratio
 * @param[in]     igrid the index of this patch in the child array
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_interpolate_face_t)(const int* mx, 
                                                          const int* my, 
                                                          const int* mz,
                                                          const int* mbc,
                                                          const int* meqn,
                                                          double qcoarse[], 
                                                          double qfine[],
                                                          const int* idir, 
                                                          const int* iface_coarse,
                                                          const int* num_neighbors,
                                                          const int* refratio, 
                                                          const int* igrid,
                                                          struct fclaw_patch_transform_data** transform_ptr);

/**
 * @brief Copies ghost data from a edge-neighboring grid on the same level
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in,out] qthis the solution of this patch
 * @param[in]     qneighbor the solution of the neighbor patch
 * @param[in]     iface the interface that the neighbor is on
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_copy_edge_t)(const int* mx, 
                                                   const int* my, 
                                                   const int* mz,
                                                   const int* mbc, 
                                                   const int* meqn,
                                                   double qthis[],
                                                   double qneighbor[], 
                                                   const int* iedge,
                                                   struct fclaw_patch_transform_data** transform_ptr);

/**
 * @brief Averages values from a edge neighboring fine grid
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in,out] qcoarse the solution of this patch
 * @param[in]     qfine the solution of the fine neighbor patch
 * @param[in]     areacoarse the area of cells in this patch
 * @param[in]     arefine the area of cells in the fine neighbor patch
 * @param[in]     manifold true if using mainifold
 * @param[in]     a_edge the edge that the fine neighbor is on
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_average_edge_t)(const int* mx, 
                                                        const int* my, 
                                                        const int* mz,
                                                        const int* mbc,
                                                        const int* meqn, 
                                                        const int* a_refratio,
                                                        double qcoarse[], 
                                                        double qfine[],
                                                        double areacoarse[], 
                                                        double areafine[],
                                                        const int* manifold,
                                                        const int* a_corner, 
                                                        struct fclaw_patch_transform_data** transform_ptr);


/**
 * @brief Copies ghost data from a corner-neighboring grid on the same level
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in,out] qthis the solution of this patch
 * @param[in]     qneighbor the solution of the neighbor patch
 * @param[in]     icorner_coarse the corner of the coarse patch to copy from
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_copy_corner_t)(const int* mx, 
                                                     const int* my, 
                                                     const int* mz,
                                                     const int* mbc,
                                                     const int* meqn, 
                                                     double qthis[],
                                                     double qneighbor[],
                                                     const int* icorner_coarse,
                                                     struct fclaw_patch_transform_data** transform_ptr);

/**
 * @brief Averages values from a corner neighboring fine grid
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in,out] qcoarse the solution of this patch
 * @param[in]     qfine the solution of the fine neighbor patch
 * @param[in]     areacoarse the area of cells in this patch
 * @param[in]     arefine the area of cells in the fine neighbor patch
 * @param[in]     manifold true if using mainifold
 * @param[in]     a_corner the corner that the neighbor is on
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_average_corner_t)(const int* mx, 
                                                        const int* my, 
                                                        const int* mz,
                                                        const int* mbc,
                                                        const int* meqn, 
                                                        const int* a_refratio,
                                                        double qcoarse[], 
                                                        double qfine[],
                                                        double areacoarse[], 
                                                        double areafine[],
                                                        const int* manifold,
                                                        const int* a_corner, 
                                                        struct fclaw_patch_transform_data** transform_ptr);

/**
 * @brief Interpolates values form a corner-neighboring coarse grid
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in]     refratio the refinement ratio
 * @param[in]     qcoarse the solution of the coarse patch
 * @param[in,out] qfine the solution of the fine patch
 * @param[in]     icorner_coarse the corner of the coarse neighbor that to interpolate from
 * @param[in]     transform_ptr Encoding for indices at block boundaries (C only).
 */
typedef void (*fclaw3d_clawpatch_fort_interpolate_corner_t)(const int* mx, 
                                                            const int* my, 
                                                            const int* mz,
                                                            const int* mbc,
                                                            const int* meqn, 
                                                            const int* refratio, 
                                                            double qcoarse[],
                                                            double qfine[], 
                                                            const int* icorner_coarse,
                                                            struct fclaw_patch_transform_data** transform_ptr);
    
/** @} */

/** @{ @name Regridding Functions *//*--------------------------------------------------*/

/**
 * @brief Tags a patch for refinement.
 * 
 * @param[in]  mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]  mbc the number of ghost cells
 * @param[in]  meqn the number of equations
 * @param[in]  xlower, ylower, zlower the coordinate of the lower bottom left corner
 * @param[in]  dx, dy, dz spacing of cells in the x, y, and z directions
 * @param[in]  blockno the block number
 * @param[in]  q the solution
 * @param[in]  tag_threshold the threshold for tagging
 * @param[in]  init_flag true if in initialization stage
 * @param[out] tag_patch true if patch should be refined
 */
typedef void (*fclaw3d_clawpatch_fort_tag4refinement_t)(const int* mx,
                                                        const int* my, 
                                                        const int* mz,
                                                        const int* mbc,
                                                        const int* meqn,
                                                        const double* xlower, 
                                                        const double* ylower, 
                                                        const double* zlower,
                                                        const double* dx, 
                                                        const double* dy, 
                                                        const double* dz,
                                                        const int* blockno,
                                                        double q[],
                                                        const double* tag_threshold,
                                                        const int* init_flag,
                                                        int* tag_patch);
/**
 * @brief Tags a quad of patches for coarsening.
 * 
 * @param[in]  mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]  mbc the number of ghost cells
 * @param[in]  meqn the number of equations
 * @param[in]  xlower, ylower, zlower the coordinate of the lower bottom left corner
 * @param[in]  dx, dy, dz spacing of cells in the x, y, and z directions
 * @param[in]  blockno the block number
 * @param[in]  q1, q2, q3, q4 the solutions on the patches
 * @param[in]  tag_threshold the threshold for tagging
 * @param[in]  init_flag true if in initialization stage
 * @param[out] tag_patch true if patches should be coarsened 
 */
typedef void (*fclaw3d_clawpatch_fort_tag4coarsening_t)(const int* mx, 
                                                        const int* my, 
                                                        const int* mz,
                                                        const int* mbc, 
                                                        const int* meqn,
                                                        double xlower[], 
                                                        double ylower[],
                                                        double zlower[],
                                                        const double* dx, 
                                                        const double* dy, 
                                                        const double* dz,
                                                        const int* blockno,
                                                        double q0[],
                                                        double q1[],
                                                        double q2[],
                                                        double q3[],
                                                        const double* tag_threshold,
                                                        const int* init_flag,
                                                        int* tag_patch);

/** 
 * @deprecated Checks if solution exceeds a threshold
 */
typedef int (*fclaw3d_clawpatch_fort_exceeds_threshold_t)(const int *blockno,
                                                          const int* meqn,
                                                          const double *qval, 
                                                          const double *qmin, 
                                                          const double *qmax,
                                                          const double quad[], 
                                                          const double *dx, 
                                                          const double *dy, 
                                                          const double *dz, 
                                                          const double *xc, 
                                                          const double *yc, 
                                                          const double *zc,
                                                          const int* ivar_variable,
                                                          const double *tag_threshold,
                                                          const int *init_flag,
                                                          const int *is_ghost);

/** 
 * @brief Averages a fine patch to a coarse patch
 * 
 * @param[in]  mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]  mbc the number of ghost cells
 * @param[in]  meqn the number of equations
 * @param[in]  qcoarse the coarse solution
 * @param[out] qfine the fine solution
 * @param[in]  areacoarse the area of the coarse cells
 * @param[in]  areafine the area of the fine cells
 * @param[in]  igrid the index of the fine patch in the siblings array
 * @param[in]  manifold true if using manifold
 */
typedef void (*fclaw3d_clawpatch_fort_interpolate2fine_t)(const int* mx, 
                                                          const int* my, 
                                                          const int* mz,
                                                          const int* mbc, 
                                                          const int* meqn,
                                                          double qcoarse[], 
                                                          double qfine[],
                                                          double areacoarse[], 
                                                          double areafine[],
                                                          const int* igrid, 
                                                          const int* manifold);
    
/** 
 * @brief Interpolates from a coarse patch to a fine patche
 * 
 * @param[in]  mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]  mbc the number of ghost cells
 * @param[in]  meqn the number of equations
 * @param[in]  qcoarse the coarse solution
 * @param[out] qfine the fine solution
 * @param[in]  areacoarse the area of the coarse cells
 * @param[in]  areafine the area of the fine cells
 * @param[in]  igrid the index of the fine patch in the siblings array
 * @param[in]  manifold true if using manifold
 */
typedef void (*fclaw3d_clawpatch_fort_average2coarse_t)(const int* mx, 
                                                        const int* my, 
                                                        const int* mz,
                                                        const int* mbc, 
                                                        const int* meqn,
                                                        double qcoarse[],
                                                        double qfine[],
                                                        double areacoarse[],
                                                        double areafine[],
                                                        const int* igrid, 
                                                        const int* manifold);
    

/** @} */

/** @{ @name time stepping *//*---------------------------------------------------------*/

/**
 * @brief Interpolates q between timesteps.
 * 
 * This only needs to interpolate the interior cells that are needed for ghost cells on
 * neighboring patches
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in]     psize    the total number cells that should be interpolated 
 * @param[in]     qcurr the current q
 * @param[in]     qlast the previous q
 * @param[out]    qinterp the inerpolated q
 * @param[in]     alpha where to interpolate between qlast and qcurr
 * @param[out]    ierror error value
 */
typedef void (*fclaw3d_clawpatch_fort_timeinterp_t)(const int *mx, 
                                                    const int* my, 
                                                    const int* mz,
                                                    const int* mbc,
                                                    const int *meqn, 
                                                    const int* psize,
                                                    double qcurr[], 
                                                    double qlast[],
                                                    double qinterp[],
                                                    const double* alpha,
                                                    const int* ierror);
    
/** @} */

/** @{ @name Parallel ghost patches *//*------------------------------------------------*/

/**
 * @brief Packs/Unpacks ghost cell data
 * 
 * @param[in]     mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of equations
 * @param[in]     mint the number of internal cells to include
 * @param[in,out] qdata the solution
 * @param[in,out] area the area of each cell
 * @param[in,out] qpack the buffer to pack from/to
 * @param[in]     psize the size of qpack buffer
 * @param[in]     packmode the packing mode (0 for packing q, 1 for unpacking q, 
 *                2 for packing q and area, and 3 for unpacking q and area) 
 * @param[out]    ierror the error value
 */
typedef void (*fclaw3d_clawpatch_fort_local_ghost_pack_t)(const int *mx, 
                                                          const int *my, 
                                                          const int *mz, 
                                                          const int *mbc,
                                                          const int *meqn, 
                                                          const int *mint,
                                                          double qdata[], 
                                                          double area[],
                                                          double qpack[], 
                                                          const int *psize,
                                                          const int *packmode, 
                                                          int *ierror);
    
/** @} */

/** @{ @name Output functions *//*------------------------------------------------------*/

/**
 * @brief Outputs the header for the time file and leaves and empty data file
 * 
 * @param[in] matname1 name of the data file
 * @param[in] matname2 name of the time file
 * @param[in] time the time
 * @param[in] meqn the number of equations
 * @param[in] maux the number of aux equations
 * @param[in] ngrids the number of grids (patches)
 */
typedef void  (*fclaw_clawpatch_fort_header_ascii_t)(const char* matname1, 
                                                     const char* matname2,
                                                     const double* time, 
                                                     const int* meqn, 
                                                     const int* maux, 
                                                     const int* ngrids);

/**
 * @brief Writes out patch data in ascii format.
 * 
 * This should append the data to the file 
 * 
 * @param[in] matname1 the name of the data file
 * @param[in] mx, my, mz the number cells in the x, y, and z directions, excluding ghost
 * @param[in] meqn the number of equations
 * @param[in] mbc the number of ghost cells
 * @param[in] xlower, ylower, zlower the coordinate of the lower bottom left corner
 * @param[in] dx, dy spacing of cells in the x, y, and z directions
 * @param[in] q the solution
 * @param[in] pach_num the patch number
 * @param[in] level the level of the patch
 * @param[in] blockno the block number
 * @param[in] mpirank the mpi rank of the patch
 */
typedef void (*fclaw3d_clawpatch_fort_output_ascii_t)(char* matname1,
                                                      int* mx,        
                                                      int* my, 
                                                      int* mz,
                                                      int* meqn,      
                                                      int* mbc,
                                                      double* xlower, 
                                                      double* ylower, 
                                                      double* zlower,
                                                      double* dx,     
                                                      double* dy, 
                                                      double* dz,
                                                      double q[],
                                                      int* patch_num, 
                                                      int* level,
                                                      int* blockno,   
                                                      int* mpirank);


/** @} */

/** @{ @name Diagnostic functions *//*--------------------------------------------------*/

/**
 * @brief Calculates the error for cells in a patch
 * 
 * @param[in] blockno the block number
 * @param[in] mx, my, mz the number cells in the x and y directions, excluding ghost
 * @param[in] mbc the number of ghost cells
 * @param[in] meqn the number of equations
 * @param[in] dx, dy, dz spacing of cells in the x and y directions
 * @param[in] xlower, ylower, zlower the coordinate of the lower left corner
 * @param[in] q the computed solution
 * @param[out] error the error for each meqn in each cell
 * @param[in] soln the exact solution
 */
typedef void (*fclaw3d_clawpatch_fort_error_t)(int* blockno, 
                                               int *mx, 
                                               int *my, 
                                               int* mz, 
                                               int *mbc,
                                               int *meqn,
                                               double *dx, 
                                               double *dy, 
                                               double* dz, 
                                               double *xlower, 
                                               double *ylower, 
                                               double *zlower,
                                               double *t, 
                                               double q[],
                                               double error[], 
                                               double soln[]);
/**
 * @brief Calculates a sum for each equation
 * 
 * @param[in] mx, my, mz the number cells in the x and y directions, excluding ghost
 * @param[in] mbc the number of ghost cells
 * @param[in] meqn the number of equations
 * @param[in] dx, dy, dz spacing of cells in the x and y directions
 * @param[in] area the area for each cell
 * @param[in] q the computed solution
 * @param[in,out] sum the current sum for each equaiton
 * @param[in,out] c_kahan the the current c values for the Kahan summation algorithm
 */
typedef void (*fclaw3d_clawpatch_fort_conscheck_t)(int *mx, 
                                                   int *my, 
                                                   int *mz, 
                                                   int* mbc, 
                                                   int* meqn,
                                                   double *dx, 
                                                   double *dy, 
                                                   double* dz,
                                                   double area[], 
                                                   double q[], 
                                                   double sum[],
                                                   double *c_kahan);

/**
 * @brief Calculates the area of a patch
 * 
 * @param[in] mx, my, mz the number cells in the x and y directions, excluding ghost
 * @param[in] mbc the number of ghost cells
 * @param[in] dx, dy, dz spacing of cells in the x and y directions
 * @param[in] area array of area values for cells
 * @return the total area of the patch
 */
typedef double (*fclaw3d_clawpatch_fort_area_t)(int *mx, 
                                                int* my, 
                                                int* mz, 
                                                int *mbc, 
                                                double *dx, 
                                                double** dy, 
                                                double *dz, 
                                                double area[]);

/**
 * @brief Calculates the error norms for a patch
 * 
 * @param[in] blockno the block number
 * @param[in] mx, my the number cells in the x and y directions, excluding ghost
 * @param[in] mbc the number of ghost cells
 * @param[in] meqn the number of equations
 * @param[in] dx, dy, dz spacing of cells in the x and y directions
 * @param[in] area array of area values for cells
 * @param[in] error error array
 * @param[out] error_norm a 2d array of  l1, l2, and inf norms for each eqn
 */
typedef void (*fclaw3d_clawpatch_fort_norm_t)(int* blockno, 
                                              int *mx, 
                                              int *my, 
                                              int* mz,
                                              int *mbc, 
                                              int *meqn,
                                              double *dx, 
                                              double *dy, 
                                              double *dz,
                                              double area[],
                                              double error[], 
                                              double error_norm[]);

/** @} */

/** 
 * @brief Checks if solution exceeds a threshold
 */
typedef int (*fclaw3d_fort_exceeds_threshold_t)(const int *blockno,
                                                const int* meqn, 
                                                const double *qval, 
                                                const double *qmin, 
                                                const double *qmax,
                                                const double quad[], 
                                                const double *dx, 
                                                const double *dy, 
                                                const double *dz,
                                                const double *xc, 
                                                const double *yc,
                                                const double *zc, 
                                                const int* ivar_threshold,
                                                const double *tag_threshold,
                                                const int    *init_flag,
                                                const int    *is_ghost);


/** @{ @name Fortran Headers *//*-------------------------------------------------------*/

/** @brief Fortran subroutine name */
#define FCLAW3D_CLAWPATCH_GET_REFINEMENT_CRITERIA \
                  FCLAW_F77_FUNC(fclaw3d_clawpatch_get_refinement_criteria, \
                                 FCLAW3D_CLAWPATCH_GET_REFINEMENT_CRITERIA)
/** @brief C declaration of fclaw3d_clawpatch_get_refinement_critera() subroutine */
int FCLAW3D_CLAWPATCH_GET_REFINEMENT_CRITERIA();


/* ------------------------------- General threshold ---------------------------------- */

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH_TAG_CRITERIA \
                  FCLAW_F77_FUNC(fclaw3d_clawpatch_tag_criteria, \
                                  FCLAW3D_CLAWPATCH_TAG_CRITERIA)

/**
 * @brief Check if the refinment threshold is exceeded
 *
 * @param[in] blockno the block number
 * @param[in] qval the 
 * @param[in] qmin the minimum q value
 * @param[in] qmax the maximum q value
 * @param[in] quad the value and adjacent values of q
 * @param[in] dx, dy the spacing in the x and y directions
 * @param[in] xc, yc the coordinate of the cell
 * @param[in] threshold the threshold
 * @param[in] init_flag true if in init stage
 * @param[in] is_ghost true if cell is a ghost cell
 * @return 1 if exceeds threshold, 0 if not, -1 if inconclusive.
 */
int FCLAW3D_CLAWPATCH_TAG_CRITERIA(const int *blockno,
                                   const double *qval, 
                                   const double *qmin, 
                                   const double *qmax,
                                   const double quad[], 
                                   const double *dx, 
                                   const double *dy, 
                                   const double *dz,
                                   const double *xc, 
                                   const double *yc,
                                   const double *zc, 
                                   const double *tag_threshold,
                                   const int *init_flag,
                                   const int *is_ghost);

/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH_EXCEEDS_THRESHOLD \
                  FCLAW_F77_FUNC(fclaw3dx_clawpatch_exceeds_threshold, \
                                  FCLAW3DX_CLAWPATCH_EXCEEDS_THRESHOLD)

/**
 * @brief Check if the refinment threshold is exceeded
 *
 * @param[in] blockno the block number
 * @param[in] qval the 
 * @param[in] qmin the minimum q value
 * @param[in] qmax the maximum q value
 * @param[in] quad the value and adjacent values of q
 * @param[in] dx, dy the spacing in the x and y directions
 * @param[in] xc, yc the coordinate of the cell
 * @param[in] threshold the threshold
 * @param[in] init_flag true if in init stage
 * @param[in] is_ghost true if cell is a ghost cell
 * @return 1 if exceeds threshold, 0 if not, -1 if inconclusive.
 */
int FCLAW3D_CLAWPATCH_EXCEEDS_THRESHOLD(const int *blockno,
                                        const int* meqn, 
                                        const double *qval, 
                                        const double *qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *dz,
                                        const double *xc, 
                                        const double *yc,
                                        const double *zc, 
                                        const int* ivar_variable,
                                        const double *tag_threshold,
                                        const int *init_flag,
                                        const int *is_ghost);


/* ----------------------------- Value threshold -------------------------------------- */
/** @brief C declaration of fclaw2d_clawpatch_value_exceeds_th() subroutine */
#define FCLAW3D_CLAWPATCH_VALUE_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw3d_clawpatch_value_exceeds_th, \
                                 FCLAW3D_CLAWPATCH_VALUE_EXCEEDS_TH)
    
/** @brief C declaration of fclaw3dx_clawpatch_value_exceeds_th() subroutine */
int FCLAW3D_CLAWPATCH_VALUE_EXCEEDS_TH(const int* blockno,
                                       const int* meqn, 
                                       const double *qval, 
                                       const double* qmin, 
                                       const double *qmax,
                                       const double quad[], 
                                       const double *dx, 
                                       const double *dy, 
                                       const double *dz,
                                       const double *xc, 
                                       const double *yc,
                                       const double *zc, 
                                       const int* ivar_variable,
                                       const double* tag_threshold,
                                       const int* init_flag,
                                       const int* is_ghost);
    
/* ----------------------------- difference threshold --------------------------------- */

/** @brief C declaration of fclaw3dx_clawpatch_difference_exceeds_th() subroutine */
#define FCLAW3D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw3d_clawpatch_difference_exceeds_th, \
                                 FCLAW3D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH)

/** @brief C declaration of fclaw3d_clawpatch_difference_exceeds_th() subroutine */
int FCLAW3D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH(const int *blockno,
                                            const int *meqn,
                                            const double *qval, 
                                            const double *qmin, 
                                            const double *qmax,
                                            const double quad[], 
                                            const double *dx, 
                                            const double *dy, 
                                            const double *dz,
                                            const double *xc, 
                                            const double *yc,
                                            const double *zc, 
                                            const int* ivar_variable, 
                                            const double *tag_threshold,
                                            const int *init_flag,
                                            const int *is_ghost);

/* --------------------------------- minmax threshold --------------------------------- */

/** @brief C declaration of fclaw3dx_clawpatch_minmax_exceeds_th() subroutine */
#define FCLAW3D_CLAWPATCH_MINMAX_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw3d_clawpatch_minmax_exceeds_th, \
                                 FCLAW3D_CLAWPATCH_MINMAX_EXCEEDS_TH)

/** @brief C declaration of fclaw3dx_clawpatch_minmax_exceeds_th() subroutine */
int FCLAW3D_CLAWPATCH_MINMAX_EXCEEDS_TH(const int *blockno,
                                        const int* meqn,
                                        const double *qval, 
                                        const double* qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *dz,
                                        const double *xc, 
                                        const double *yc,
                                        const double *zc, 
                                        const int* ivar_variable,
                                        const double *tag_threshold,                        
                                        const int *init_flag,
                                        const int *is_ghost);

/* ------------------------------- gradient threshold --------------------------------- */
/** @brief C declaration of fclaw3dx_clawpatch_gradient_exceeds_th() subroutine */
#define FCLAW3D_CLAWPATCH_GRADIENT_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw3d_clawpatch_gradient_exceeds_th, \
                                 FCLAW3D_CLAWPATCH_GRADIENT_EXCEEDS_TH)

/** @brief C declaration of fclaw3dx_clawpatch_gradient_exceeds_th() subroutine */
int FCLAW3D_CLAWPATCH_GRADIENT_EXCEEDS_TH(const int *blockno,
                                          const int* meqn,
                                          const double *qval, 
                                          const double* qmin, 
                                          const double *qmax,
                                          const double quad[], 
                                          const double *dx, 
                                          const double *dy, 
                                          const double *dz,
                                          const double *xc, 
                                          const double *yc,
                                          const double *zc, 
                                          const int* ivar_variable,
                                          const double *tag_threshold,
                                          const int *init_flag,
                                          const int *is_ghost);
    

/* ------------------------------- user threshold --------------------------------- */
/** Fortran subroutine name */
#define FCLAW3D_USER_EXCEEDS_TH FCLAW_F77_FUNC( \
                             fclaw3d_user_exceeds_th, \
                             FCLAW3D_USER_EXCEEDS_TH)

/** @brief C declaration of user_exceeds_th() subroutine */
int FCLAW3D_USER_EXCEEDS_TH(const int *blockno,
                            const int* meqn,
                            const double *qval, 
                            const double* qmin, 
                            const double *qmax,
                            const double quad[], 
                            const double *dx, 
                            const double *dy, 
                            const double *dz,
                            const double *xc, 
                            const double *yc,
                            const double *zc, 
                            const int* ivar_variable,
                            const double *tag_threshold,
                            const int *init_flag,
                            const int *is_ghost);

/* -------------------------- User convenience headers -------------------------------- */

/** Fortran subroutine name */
#define FCLAW3D_USER_TAG4REFINEMENT FCLAW_F77_FUNC( \
                 fclaw3d_user_tag4refinement, FCLAW3D_USER_TAG4REFINEMENT)
/** 
 * @brief C declaration of user defined tag4refinement subroutine, see 
 * ::clawpatch_fort_tag4refinement_t
 */
void FCLAW3D_USER_TAG4REFINEMENT(const int* mx,
                                 const int* my,
                                 const int* mz,
                                 const int* mbc,
                                 const int* meqn,
                                 const double* xlower, 
                                 const double* ylower,
                                 const double* zlower,
                                 const double* dx, 
                                 const double* dy,
                                 const double* dz,
                                 const int* blockno,
                                 double q[],
                                 const double* tag_threshold,
                                 const int* init_flag,
                                 int* tag_patch);

/** Fortran subroutine name */
#define FCLAW3D_USER_TAG4COARSENING FCLAW_F77_FUNC( \
            fclaw3d_user_tag4coarsening, FCLAW3D_USER_TAG4COARSENING)
/** 
 * @brief C declaration of user defined tag4coarsening subroutine, see 
 * ::clawpatch_fort_tag4coarsening_t
 */
void FCLAW3D_USER_TAG4COARSENING(const int* mx, 
                                 const int* my,
                                 const int* mz,
                                 const int* mbc, 
                                 const int* meqn,
                                 const double* xlower, 
                                 const double* ylower,
                                 const double* zlower,
                                 const double* dx, 
                                 const double* dy,
                                 const double* dz,
                                 const int* blockno,
                                 double q0[],
                                 double q1[],
                                 double q2[],
                                 double q3[],
                                 const double* tag_threshold,
                                 const int* initflag,
                                 int* tag_patch);
    

/* ----------------------------- interpolation/coarsening ----------------------------- */

/** Fortran subroutine name */
#define FCLAW3D_USER_INTERPOLATE2FINE FCLAW_F77_FUNC(\
                  fclaw3d_user_interpolate2fine, FCLAW3D_USER_INTERPOLATE2FINE)
/**
 * @brief C declaration of user defined interpolate2fine subroutine,
 * see ::clawpatch_fort_interpolate2fine_t
 */
void FCLAW3D_USER_INTERPOLATE2FINE(const int* mx,
                                   const int* my,
                                   const int* mz,
                                   const int* mbc,
                                   const int* meqn, 
                                   double qcoarse[], 
                                   double qfine[],
                                   double areacoarse[], 
                                   double areafine[],
                                   const int* igrid,
                                   const int* manifold);

/** Fortran subroutine name */
#define FCLAW3D_USER_AVERAGE2COARSE FCLAW_F77_FUNC(\
             fclaw3d_user_average2coarse, FCLAW3D_USER_AVERAGE2COARSE)
/**
 * @brief C declaration of user defined average2coarse subroutine,
 * see ::clawpatch_fort_average2coarse_t
 */
void FCLAW3D_USER_AVERAGE2COARSE(const int* mx,
                                 const int* my,
                                 const int* mz,
                                 const int* mbc,
                                 const int* meqn,
                                 double qcoarse[],
                                 double qfine[],
                                 double areacoarse[],
                                 double areafine[],
                                 const int* igrid, 
                                 const int* manifold);
    
/** @} */

#ifdef __cplusplus
}
#endif

#endif
