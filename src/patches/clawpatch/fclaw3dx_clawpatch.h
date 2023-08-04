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

#ifndef FCLAW3DX_CLAWPATCH_H
#define FCLAW3DX_CLAWPATCH_H

#include <fclaw2d_defs.h>      /* Needed to get correction def. of PATCHDIM */

#include <forestclaw2d.h>       /* Need patch callback def */

#include <fclaw_clawpatch_enums.h>

#include <fclaw2d_clawpatch.h>

#include <fclaw3dx_clawpatch_fort.h>
#include <fclaw_clawpatch_diagnostics.h>


#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

/** 
 *  @file
 *  Clawpatch related structures and routines
 */


/* ---------------------------- Virtual table ------------------------------------ */
/* members of this structure provide the only access to above functions */

/**
 * @brief Initialize the clawpatch vtable global variable
 * 
 * @param claw_version the version of clawpack (4 for 4.6, 5 for 5)
 */
void fclaw3dx_clawpatch_vtable_initialize(struct fclaw2d_global *glob, int claw_version);

/* -------------------------------- time stepping ----------------------------------- */

/* Called in step2 (clawpack 4.6 and clawpack 5.0) */

/**
 * @brief Save the current timestep
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 */
void fclaw3dx_clawpatch_save_current_step(struct fclaw2d_global* glob,
                                          struct fclaw2d_patch* patch);


/* ------------------------------- Misc access functions ------------------------------ */

/**
 * @brief Get the grid data for a specific patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] mx, my, mz the number of cells in the x, y, and z directions
 * @param[out] mbc the number of ghost cells
 * @param[out] xlower, ylower, zlower the lower bottom left coordinate of the patch
 * @param[out] dx, dy, dz the spacings in the x, y, and z directions
 */
void fclaw3dx_clawpatch_grid_data(struct fclaw2d_global* glob,
                                  struct fclaw2d_patch* patch,
                                  int* mx, 
                                  int* my, 
                                  int* mz, 
                                  int* mbc,
                                  double* xlower, 
                                  double* ylower,
                                  double* zlower, 
                                  double* dx, 
                                  double* dy, 
                                  double* dz);


void fclaw3d_clawpatch_grid_data(struct fclaw2d_global* glob,
                                  struct fclaw2d_patch* patch,
                                  int* mx, 
                                  int* my, 
                                  int* mz, 
                                  int* mbc,
                                  double* xlower, 
                                  double* ylower,
                                  double* zlower, 
                                  double* dx, 
                                  double* dy, 
                                  double* dz);


/**
 * @brief Get the area for each cell of the patch
 * 
 * @param this_patch the patch to get the are for
 * @return double* the area array
 */
double* fclaw3d_clawpatch_get_volume(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* this_patch);

/**
 * @brief Get the scalar metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] area the area of each cell in the patch
 * @param[out] edgelengths the edge lenghts for each cell. 
 *             An array of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the the bottom edge length for cell (i,j).
 *             (i,j,2) contains the the left edge length for cell (i,j).
 * @param[out] curvature the curvature for each cell in the patch
 */
void fclaw3d_clawpatch_metric_scalar(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch,
                                     double **volume,
                                     double **faceareas);

/**
 * @brief Get the vector metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xnormals, ynormals the normals of each face of the patch
 *             Arrays of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the normal for the bottom edge of cell (i,j).
 *             (i,j,2) contains the normal for the left edge of cell (i,j).
 * @param[out] xtangents, ytangents the tangents of each face of the patch
 *             Arrays of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the tangent for the bottom edge of cell (i,j).
 *             (i,j,2) contains the tangent for the left edge of cell (i,j).
 * @param[out] surfnormals the surface normal for each cell center of the patch
 *             An array of dimension(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)
 */
void fclaw3d_clawpatch_metric_basis(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch* patch,
                                      double **xrot, 
                                      double **yrot,
                                      double **zrot);

/**
 * @brief Get the mesh metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xp, yp, zp the coordinates of the cell centers
 * @param[out] xd, yd, zd the coordinates of the nodes
 * @param[out] area the area of each cell
 */
void fclaw3d_clawpatch_mesh_data(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* patch,
                                 double **xp, 
                                 double **yp, 
                                 double **zp,
                                 double **xd, 
                                 double **yd, 
                                 double **zd,
                                 double **volume,
                                 double **faceareas);

/**
 * @brief Get the solution data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @param[out] q the solution array
 * @param[out] meqn the number of equations
 */
void fclaw3dx_clawpatch_soln_data(struct fclaw2d_global* glob,
                                  struct fclaw2d_patch* patch,
                                  double **q, 
                                  int* meqn);
/**
 * @brief Get the aux data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @param[out] aux the aux array
 * @param[out] maux the number of equations
 */
void fclaw3dx_clawpatch_aux_data(struct fclaw2d_global *glob,
                                 struct fclaw2d_patch *patch,
                                 double **aux, 
                                 int* maux);

/**
 * @brief Get the rhs data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @param[out] rhs the rhs array
 * @param[out] mfields the number fields
 */
void fclaw3dx_clawpatch_rhs_data(struct fclaw2d_global* glob,
                                 fclaw2d_patch_t* patch,
                                 double **rhs, 
                                 int *mfields);
/**
 * @brief Get the error data for elliptic problems
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @param[out] err the error array
 * @param[out] mfields the number fields
 */
void fclaw3dx_clawpatch_elliptic_error_data(struct fclaw2d_global* glob,
                                            struct fclaw2d_patch* patch,
                                            double **err, 
                                            int *mfields);
/**
 * @brief Get the solution data for elliptic problems
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @param[out] soln the solution array
 * @param[out] mfields the number fields
 */
void fclaw3dx_clawpatch_elliptic_soln_data(struct fclaw2d_global* glob,
                                           struct fclaw2d_patch* patch,
                                           double **soln, 
                                           int *mfields);

/**
 * @brief Get the solution data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @return double* the solution array
 */
double* fclaw3dx_clawpatch_get_q(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* patch);

/**
 * @brief Get the error data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @return double* the error array
 */
double* fclaw3dx_clawpatch_get_error(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch);

/**
 * @brief Get the exact solution data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @return double* the exact solution array
 */
double* fclaw3dx_clawpatch_get_exactsoln(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* patch);

/**
 * @brief Get the size of a solution array for a patch
 * 
 * @param glob the global context
 * @return size_t the size
 */
size_t fclaw3dx_clawpatch_size(struct fclaw2d_global *glob);


/**
 * @brief Get the user data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return void* the pointer
 */
void* fclaw3dx_clawpatch_get_user_data(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch* patch);

/**
 * @brief Set the user data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @param udata the user data pointer
 */
void fclaw3dx_clawpatch_set_user_data(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch* patch,
                                      void* udata);
/**
 * @brief Get the solver data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return void* the pointer
 */
void* fclaw3dx_clawpatch_get_solver_data(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* patch);

/**
 * @brief Set the solver data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @param sdata the solver data pointer
 */
void fclaw3dx_clawpatch_set_solver_data(struct fclaw2d_global* glob,
                                        struct fclaw2d_patch* patch,
                                        void* sdata);


/* These should be renamed to time_interp data */

/**
 * @brief Get an interpolated solution
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] time_interp true if interpolated grid data should be returned
 * @param[out] q the interpolated solution
 * @param[out] meqn the number of equations
 */
void fclaw3dx_clawpatch_timesync_data(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch* patch,
                                      int time_interp,
                                      double **q, 
                                      int* meqn);

/**
 * @brief Get an interpolated solution
 * 
 * @param glob the global context
 * @param this_patch the patch context
 * @param time_interp true if interpolated grid data should be returned
 * @return the interpolated solution
 */
double* fclaw3dx_clawpatch_get_q_timesync(struct fclaw2d_global* glob,
                                          struct fclaw2d_patch* patch,
                                          int time_interp);
/**
 * @brief Get the registers for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return the registers
 */
struct fclaw3dx_clawpatch_registers* 
fclaw3dx_clawpatch_get_registers(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* patch);

#ifdef __cplusplus
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_H */
