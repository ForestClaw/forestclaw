/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW2D_CLAWPATCH_H
#define FCLAW2D_CLAWPATCH_H

#include <fclaw2d_defs.h>      /* Needed to get correction def. of PATCHDIM */

#include <forestclaw2d.h>       /* Need patch callback def */

#include <fclaw2d_clawpatch_fort.h>
#include <fclaw2d_clawpatch_conservation.h>
#include <fclaw2d_clawpatch_diagnostics.h>


#ifdef __cplusplus
extern "C"
{
#endif

/** 
 *  @file
 *  Clawpatch related structures and routines
 */

#if 0
/* Fix syntax highlighting */
#endif


/**
 * @brief fclaw2d_clawpatch_vtable type
 */
typedef struct fclaw2d_clawpatch_vtable fclaw2d_clawpatch_vtable_t;

/* --------------------------------- Typedefs ----------------------------------------- */
/**
 * @brief Sets a pointer to user data for a specific patch
 * 
 * @param[in] glob the global context
 * @param[in] patch the patch context
 * @param[in] user the pointer to the user data
 */
typedef void (*clawpatch_set_user_data_t)(struct fclaw2d_global *glob, 
                                          struct fclaw2d_patch *patch,
                                          void* user);

/**
 * @brief Packs/Unpacks the fclaw2d_clawpatch_registers struct for ghost patches
 * 
 * @param[in] glob the global context
 * @param[in] patch the patch context
 * @param[in,out] qpack the buffer
 * @param[in] frsize the size of the buffer
 * @param[in] packmode enum fclaw2d_clawpatch_packmode enum
 * @param[out] ierror error value
 */
typedef void (*clawpatch_time_sync_pack_registers_t)(struct fclaw2d_global *glob,
                                                     struct fclaw2d_patch *this_patch,
                                                     double *qpack,
                                                     int frsize, 
                                                     fclaw2d_clawpatch_packmode_t packmode,
                                                     int *ierror);


/**
 * @brief Packs/Unpacks ghost cell data for an aux array
 * 
 * @param[in]     glob the global context
 * @param[in]     patch the patch context
 * @param[in]     mint the number of internal cells to pack
 * @param[in,out] qpack the the buffer
 * @param[in]     extrasize the size of the buffer
 * @param[in]     packmode the packing mode (0 for packing aux, 1 for unpacking aux)
 * @param[out]    ierror error value
 */
typedef void (*clawpatch_local_ghost_pack_aux_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *patch,
                                                 const int mint,
                                                 double qpack[], int extrasize,
                                                 int packmode, int* ierror);
    
/**
 * @brief Packs/Unpacks the fclaw2d_clawpatch_registers struct for ghost patches
 * 
 * @param[in]     glob the global context
 * @param[in]     patch the patch context
 * @param[in,out] qpack the buffer
 * @param[in]     the size of the buffer
 * @param[out]    ierror error value
 */
typedef void (*clawpatch_fort_local_ghost_pack_registers_t)(struct fclaw2d_global *glob,
                                                            struct fclaw2d_patch *patch,
                                                            double qpack[], int frsize,
                                                            int* ierror);
/* ------------------------------ typedefs - output ----------------------------------- */

/**
 * @brief Outputs a time header
 * 
 * @param[in] glob the global context
 * @param[in] iframe the frame
 */
typedef void (*clawpatch_time_header_t)(struct fclaw2d_global* glob, int iframe);


/* ---------------------------- typedefs - diagnostics -------------------------------- */

/**
 * @brief Fills in a user defined error_data structure for a specific patch.
 * 
 * @param[in] glob glob the global context
 * @param[in] patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in,out] error_data a user defined structure
 * 
 */
typedef void (*clawpatch_diagnostics_cons_t)(struct fclaw2d_global *glob,
                                             struct fclaw2d_patch *patch,
                                             int blockno,
                                             int patchno,
                                             void *error_data);

/**
 * @brief Fills in a user defined error_data structure for a specific patch.
 * 
 * @param[in] glob glob the global context
 * @param[in] patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in,out] error_data a user defined structure
 * 
 */
typedef void (*clawpatch_diagnostics_error_t)(struct fclaw2d_global *glob,
                                              struct fclaw2d_patch *patch,
                                              int blockno,
                                              int patchno,
                                              void *error_data);


/* ---------------------------- Virtual table ------------------------------------ */
/* members of this structure provide the only access to above functions */

/**
 * @brief Initialize the clawpatch vtable global variable
 * 
 * @param claw_version the version of clawpack (4 for 4.6, 5 for 5)
 */
void fclaw2d_clawpatch_vtable_initialize(int claw_version);

/**
 * @brief Get a pointer to a clawpatch vtable global variable
 * 
 * @return fclaw2d_clawpatch_vtable_t* the vtable
 */
fclaw2d_clawpatch_vtable_t* fclaw2d_clawpatch_vt();

/**
 * @brief vtable for clawpatch related functions
 */
struct fclaw2d_clawpatch_vtable
{
    /**
     * @brief Function that allows the user to set a user data pointer
     */
    clawpatch_set_user_data_t              set_user_data;

    /** @{ @name Ghost Filling Functions */

    /** Copies ghost data from a face neighboring grid on the same level */
    clawpatch_fort_copy_face_t             fort_copy_face;
    /** Averages values from a face neighboring fine grid */
    clawpatch_fort_average_face_t          fort_average_face;
    /** Interpolates values form a face neighboring coarse grid */
    clawpatch_fort_interpolate_face_t      fort_interpolate_face;
    /** Copies ghost data from a corner neighboring grid on the same level */
    clawpatch_fort_copy_corner_t           fort_copy_corner;
    /** Averages values from a corner neighboring fine grid */
    clawpatch_fort_average_corner_t        fort_average_corner;
    /** Interpolates values form a corner neighboring coarse grid */
    clawpatch_fort_interpolate_corner_t    fort_interpolate_corner;

    /** @} */

    /** @{ @name Regridding Functions */

    /** Tags a patch for refinement. */
    clawpatch_fort_tag4refinement_t        fort_tag4refinement;
    /** Tags a quad of patches for coarsening. */
    clawpatch_fort_tag4coarsening_t        fort_tag4coarsening;
    /** @deprecated Checks if solution exceeds a threshold */
    clawpatch_fort_exceeds_threshold_t     fort_user_exceeds_threshold;

    /** Averages a fine patches to a coarse patch */
    clawpatch_fort_average2coarse_t        fort_average2coarse;
    /** Interpolates from a coarse patch to a fine patche */
    clawpatch_fort_interpolate2fine_t      fort_interpolate2fine;

    /** @} */

    /** @{ @name Conservation Update */

    /** Adds fine grid corrections to coarse grid. */
    clawpatch_fort_time_sync_f2c_t         fort_time_sync_f2c;
    /** Adds wave corrections at same level interfaces. */
    clawpatch_fort_time_sync_samesize_t    fort_time_sync_samesize;
    /** Packs/Unpacks the fclaw2d_clawpatch_registers struct for ghost patches */
    clawpatch_time_sync_pack_registers_t   time_sync_pack_registers;

    /** @} */

    /** @{ @name Output Functions (ascii) */

    /** Outputs a time header */
    clawpatch_time_header_t                time_header_ascii;
    /** Outputs a time header */
    clawpatch_fort_header_ascii_t          fort_header_ascii;

    /** Called for every patch when outputing ascii */
    fclaw2d_patch_callback_t               cb_output_ascii;    
    /** Outputs patch data in ascii */
    clawpatch_fort_output_ascii_t          fort_output_ascii;

    /** @} */

    /** @{ @name Time interpolation functions */

    /** Interpolates q between timesteps */
    clawpatch_fort_timeinterp_t            fort_timeinterp;

    /** @} */

    /** @{ @name Ghost Patch Functions */
    
    /** Packs/Unpacks ghost cell data */
    clawpatch_fort_local_ghost_pack_t      fort_local_ghost_pack;
    /** Packs/Unpacks ghost cell data for an aux array */
    clawpatch_local_ghost_pack_aux_t       local_ghost_pack_aux;

    /** @} */

    /** @{ @name Diagnostic Functions */

    /** Fills in a user defined error_data structure for a specific patch. */
    clawpatch_diagnostics_cons_t           conservation_check;
    /** Fills in a user defined error_data structure for a specific patch. */
    clawpatch_diagnostics_error_t          compute_error;

    /** Calculates the error for cells in a patch */
    clawpatch_fort_error_t                 fort_compute_patch_error;
    /** Calculates a sum for each equation */
    clawpatch_fort_conscheck_t             fort_conservation_check;
    /** Calculates the error norms for a patch */
    clawpatch_fort_norm_t                  fort_compute_error_norm;
    /** Calculates the area of a patch */
    clawpatch_fort_area_t                  fort_compute_patch_area;

    /** @} */

    /** @{ @name Diagnostics */

    /** Whether or not this vtable is set */
    int is_set; 

    /** @} */
};


/* -------------------------------- time stepping ----------------------------------- */

/* Called in step2 (clawpack 4.6 and clawpack 5.0) */

/**
 * @brief Save the current timestep
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 */
void fclaw2d_clawpatch_save_current_step(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* this_patch);


/* ------------------------------- Misc access functions ------------------------------ */

/**
 * @brief Get the grid data for a specific patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] mx, my the number of cells in the x and y directions
 * @param[out] mbc the number of ghost cells
 * @param[out] xlower, ylower the lower left coordinate of the patch
 * @param[out] dx, dy the spacings in the x and y directions
 */
void fclaw2d_clawpatch_grid_data(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* patch,
                                 int* mx, 
                                 int* my, 
                                 int* mbc,
                                 double* xlower, 
                                 double* ylower,
                                 double* dx, 
                                 double* dy);
#if 0
void fclaw2d_clawpatch_grid_data(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* patch,
                                 int* mx, int* my, int* mz, 
                                 int* mbc,
                                 double* xlower, double* ylower,
                                 double* zlower, 
                                 double* dx, double* dy, double* dz);
#endif

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
void fclaw2d_clawpatch_metric_scalar(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch,
                                     double **area,
                                     double **edgelengths,
                                     double **curvature);

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
void fclaw2d_clawpatch_metric_vector(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch,
                                     double **xnormals, 
                                     double **ynormals,
                                     double **xtangents, 
                                     double **ytangents,
                                     double **surfnormals);

/**
 * @brief Get the mesh metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xp, yp, zp the coordinates of the cell centers
 * @param[out] xd, yd, zd the coordinates of the nodes
 * @param[out] area the area of each cell
 */
void fclaw2d_clawpatch_metric_data(struct fclaw2d_global* glob,
                                   struct fclaw2d_patch* patch,
                                   double **xp, 
                                   double **yp, 
                                   double **zp,
                                   double **xd, 
                                   double **yd, 
                                   double **zd,
                                   double **area);
/**
 * @brief Get the mesh metrics of a patch
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
 * @param[out] edgelengths the edge lenghts for each cell. 
 *             An array of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the the bottom edge length for cell (i,j).
 *             (i,j,2) contains the the left edge length for cell (i,j).
 * @param[out] curvature the curvature for each cell in the patch
 */
void fclaw2d_clawpatch_metric_data2(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* patch,
                                    double **xnormals, 
                                    double **ynormals,
                                    double **xtangents, 
                                    double **ytangents,
                                    double **surfnormals, 
                                    double ** edgelengths,
                                    double **curvature);

/**
 * @brief Get the area for each cell of the patch
 * 
 * @param this_patch the patch to get the are for
 * @return double* the area array
 */
double* fclaw2d_clawpatch_get_area(struct fclaw2d_global* glob,
                                   struct fclaw2d_patch* patch);
/**
 * @brief Get the solution data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @param[out] q the solution array
 * @param[out] meqn the number of equations
 */
void fclaw2d_clawpatch_soln_data(struct fclaw2d_global* glob,
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
void fclaw2d_clawpatch_aux_data(struct fclaw2d_global *glob,
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
void fclaw2d_clawpatch_rhs_data(struct fclaw2d_global* glob,
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
void fclaw2d_clawpatch_elliptic_error_data(struct fclaw2d_global* glob,
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
void fclaw2d_clawpatch_elliptic_soln_data(struct fclaw2d_global* glob,
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
double* fclaw2d_clawpatch_get_q(struct fclaw2d_global* glob,
                                struct fclaw2d_patch* patch);

/**
 * @brief Get the error data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @return double* the error array
 */
double* fclaw2d_clawpatch_get_error(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* patch);

/**
 * @brief Get the exact solution data for a patch
 * 
 * @param[in]  glob the global context
 * @param[in]  this_patch the patch context
 * @return double* the exact solution array
 */
double* fclaw2d_clawpatch_get_exactsoln(struct fclaw2d_global* glob,
                                        struct fclaw2d_patch* patch);

/**
 * @brief Get the size of a solution array for a patch
 * 
 * @param glob the global context
 * @return size_t the size
 */
size_t fclaw2d_clawpatch_size(struct fclaw2d_global *glob);


/**
 * @brief Get the user data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return void* the pointer
 */
void* fclaw2d_clawpatch_get_user_data(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch* patch);

/**
 * @brief Set the user data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @param udata the user data pointer
 */
void fclaw2d_clawpatch_set_user_data(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch,
                                     void* udata);
/**
 * @brief Get the solver data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return void* the pointer
 */
void* fclaw2d_clawpatch_get_solver_data(struct fclaw2d_global* glob,
                                        struct fclaw2d_patch* patch);

/**
 * @brief Set the solver data pointer for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @param sdata the solver data pointer
 */
void fclaw2d_clawpatch_set_solver_data(struct fclaw2d_global* glob,
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
void fclaw2d_clawpatch_timesync_data(struct fclaw2d_global* glob,
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
double* fclaw2d_clawpatch_get_q_timesync(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* patch,
                                         int time_interp);
/**
 * @brief Get the registers for a patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return the registers
 */
struct fclaw2d_clawpatch_registers* 
fclaw2d_clawpatch_get_registers(struct fclaw2d_global* glob,
                                struct fclaw2d_patch* patch);

#ifdef __cplusplus
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_H */
