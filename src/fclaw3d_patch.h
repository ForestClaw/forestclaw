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
/** 
 * @file
 * @brief Patch related functions and typedefs
 */


#ifndef FCLAW3D_PATCH_H
#define FCLAW3D_PATCH_H

#include <forestclaw3d.h>  /* Contains definition of patch-iterator callback */
#include <fclaw2d_patch.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


struct fclaw_global;
struct fclaw_domain;
struct fclaw3d_patch;


/* ------------------------------------------------------------------------------------ */
///                         @name Creating/Deleting Patches
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * DEPRECATED
 * @deprecated NOT USED
 */
void fclaw3d_patch_reset_data(struct fclaw_global* glob,
                              struct fclaw_patch* old_patch,
                              struct fclaw_patch* new_patch,
                              int blockno,int old_patchno, int new_patchno);


/**
 * @brief Deallocate the user data pointer for a patch
 * 
 * @param[in] glob the global context 
 * @param[in,out] patch the patch context, user data pointer is set to NULL on return
 */
void fclaw3d_patch_data_delete(struct fclaw_global *glob,
                               struct fclaw_patch *patch);

/**
 * @brief Construct a new patch object
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in,out] user user data pointer
 */
void fclaw3d_patch_build(struct fclaw_global *glob,
                         struct fclaw_patch *this_patch,
                         int blockno,
                         int patchno,
                         void *user);

/**
 * @brief Construct a new patch object from a set of fine patches
 * 
 * @param[in] glob the global context
 * @param[in] fine_patches the fine patch contexts
 * @param[in,out] coarse_patches the coarse patch context
 * @param[in] blockno the block number
 * @param[in] coarse_patchno coarse patch number
 * @param[in] fine0_patchno first fine patch number
 * @param[in] build_mode the build mode
 */
void fclaw3d_patch_build_from_fine(struct fclaw_global *glob,
                                   struct fclaw_patch *fine_patches,
                                   struct fclaw_patch *coarse_patch,
                                   int blockno,
                                   int coarse_patchno,
                                   int fine0_patchno,
                                   fclaw_build_mode_t build_mode);


///@}
/* ------------------------------------------------------------------------------------ */
///                         @name Solver Specific Functions
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Initialize patch data for a solver
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 */
void fclaw3d_patch_initialize(struct fclaw_global *glob,
                              struct fclaw_patch *this_patch,
                              int blockno,
                              int patchno);

/**
 * @brief Initialize boundary conditions
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in] t the time
 * @param[in] dt the timestep
 * @param[in] intersects_bc array of values for each face, true if physical boundary
 * @param[in] time_interp true if in time interpolation stage (not global)
 */
void fclaw3d_patch_physical_bc(struct fclaw_global *glob,
                               struct fclaw_patch *this_patch,
                               int blockno,
                               int patchno,
                               double t,
                               double dt,
                               int *intersects_bc,
                               int time_interp);

/**
 * @brief Advance a patch with a single time step
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number 
 * @param[in] patchno the patch number
 * @param[in] t the time
 * @param[in] dt the timestep
 * @param[in] user pointer to the ::fclaw3d_single_step_buffer_data struct (used in cudaclaw)
 * @return double the maxcfl
 */
double fclaw3d_patch_single_step_update(struct fclaw_global *glob,
                                        struct fclaw_patch *this_patch,
                                        int blockno,
                                        int patchno,
                                        double t,
                                        double dt, void* user);

/**
 * @brief Set the right hand side for a patch
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number 
 * @param[in] patchno the patch number
 */
void fclaw3d_patch_set_rhs(struct fclaw_global *glob,
                           struct fclaw_patch *patch,
                           int blockno,
                           int patchno);


///@}
/* ------------------------------------------------------------------------------------ */
///                              @name Time Stepping
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Restores a previously saved solution
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 */
void fclaw3d_patch_restore_step(struct fclaw_global* glob,
                                struct fclaw_patch* this_patch);

/**
 * @brief Saves the current solution for later use
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 */
void fclaw3d_patch_save_step(struct fclaw_global* glob,
                             struct fclaw_patch* this_patch);

/**
 * @brief Sets up interpolated values for a patch 
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] alpha the alpha value, with 0 being the last time step and 1 being the current time step
 */
void fclaw3d_patch_setup_timeinterp(struct fclaw_global *glob,
                                    struct fclaw_patch *this_patch,
                                    double alpha);

///@}
/* ------------------------------------------------------------------------------------ */
///                              @name Ghost Filling
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Copies ghost data from a face-neighboring grid on the same level
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch this patch context
 * @param[in] neighbor_patch the neighbor patch context
 * @param[in] iface the interface that the neighbor patch is on
 * @param[in] time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in] tranform_data the tranform data for the neighbor's coordinate system
 */
void fclaw3d_patch_copy_face(struct fclaw_global* glob,
                             struct fclaw_patch *this_patch,
                             struct fclaw_patch *neighbor_patch,
                             int iface,
                             int time_interp,
                             struct fclaw_patch_transform_data *transform_data);

/**
 * @brief Averages values from a face-neighboring fine grid
 * 
 * @param[in]     glob the global context
 * @param[in,out] coarse_patch this patch context
 * @param[in]     fine_patch the fine patch context
 * @param[in]     idir Face orientation - 0 for x-faces; 1 for y-faces [0-1]
 * @param[in]     ifaced_coarse the interface of the fine neighbor patch
 * @param[in]     num_neighbors the number of neighbors
 * @param[in]     refine_factor the refinement factor (number of neighbors)
 * @param[in]     refratio the refinement ratio
 * @param[in]     time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in]     igrid the index of the fine neighbor in the child array
 * @param[in]     tranform_data the tranform data for the neighbor's coordinate system
 */
void fclaw3d_patch_average_face(struct fclaw_global* glob,
                                struct fclaw_patch *coarse_patch,
                                struct fclaw_patch *fine_patch,
                                int idir,
                                int iface_coarse,
                                int refine_factor,
                                int refratio,
                                int time_interp,
                                int igrid,
                                struct fclaw_patch_transform_data* transform_data);

/**
 * @brief Interpolates values from a face-neighboring coarse grid
 * 
 * @param[in]     glob the global context
 * @param[in]     coarse_patch this patch context
 * @param[in,out] fine_patch the fine patch context
 * @param[in]     idir Face orientation - 0 for x-faces; 1 for y-faces [0-1]
 * @param[in]     ifaced_coarse the interface of the fine neighbor patch
 * @param[in]     num_neighbors the number of neighbors
 * @param[in]     refine_factor the refinement factor (number of neighbors)
 * @param[in]     refratio the refinement ratio
 * @param[in]     time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in]     igrid the index of the fine neighbor in the child array
 * @param[in]     tranform_data the tranform data for the neighbor's coordinate system
 */
void fclaw3d_patch_interpolate_face(struct fclaw_global* glob,
                                    struct fclaw_patch *coarse_patch,
                                    struct fclaw_patch *fine_patch,
                                    int idir,
                                    int iside,
                                    int RefineFactor,
                                    int refratio,
                                    int time_interp,
                                    int igrid,
                                    struct fclaw_patch_transform_data* transform_data);

/* Do we want fclaw3d_patch_copy/average/interpolate_edge as well?
 * Feel free to add. */

/**
 * @brief Copies values from a corner-neighboring grid
 * 
 * @param[in]     glob the global context
 * @param[in,out] this_patch this patch context
 * @param[in]     neighbor_patch the neighbor patch context
 * @param[in]     this_blockno the block number of this patch
 * @param[in]     neighbor_blockno the block number of the neighbor patch
 * @param[in]     is_block_corner true if corner is on the corner of a block
 * @param[in]     icorner the corner that the neighboring patch is on
 * @param[in]     time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in]     tranform_data the tranform data for the neighbor's coordinate system
 */
void fclaw3d_patch_copy_corner(struct fclaw_global* glob,
                               struct fclaw_patch *this_patch,
                               struct fclaw_patch *neighbor_patch,
                               int this_blockno,
                               int neighbor_blockno,
                               int is_block_corner,
                               int icorner,
                               int time_interp,
                               struct fclaw_patch_transform_data *transform_data);

/**
 * @brief Averages values from a corner-neighboring fine grid
 * 
 * @param[in]     glob the global context
 * @param[in,out] coarse_patch the coarse patch context
 * @param[in]     fine_patch the fine patch context
 * @param[in]     coarse_blockno the block number of the coarse patch
 * @param[in]     fine_blockno the block number of the fine patch
 * @param[in]     is_block_corner true if corner is on the corner of a block
 * @param[in]     icorner the corner of the coarse patch that the fine patch is on
 * @param[in]     time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in]     tranform_data the tranform data for the neighbor's coordinate system
 */
void fclaw3d_patch_average_corner(struct fclaw_global* glob,
                                  struct fclaw_patch *coarse_patch,
                                  struct fclaw_patch *fine_patch,
                                  int coarse_blockno,
                                  int fine_blockno,
                                  int is_block_corner,
                                  int coarse_corner,
                                  int time_interp,
                                  struct fclaw_patch_transform_data* transform_data);

/**
 * @brief Interpolates values from a corner-neighboring coarse grid
 * 
 * @param[in]     glob the global context
 * @param[in]     coarse_patch the coarse patch context
 * @param[in,out] fine_patch the fine patch context
 * @param[in]     coarse_blockno the block number of the coarse patch
 * @param[in]     fine_blockno the block number of the fine patch
 * @param[in]     is_block_corner true if corner is on the corner of a block
 * @param[in]     icorner the corner of the coarse patch that the fine patch is on
 * @param[in]     time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in]     tranform_data the tranform data for the neighbor's coordinate system
 */
void fclaw3d_patch_interpolate_corner(struct fclaw_global* glob,
                                      struct fclaw_patch* coarse_patch,
                                      struct fclaw_patch* fine_patch,
                                      int coarse_blockno,
                                      int fine_blockno,
                                      int is_block_corner,
                                      int coarse_corner,
                                      int time_interp,
                                      struct fclaw_patch_transform_data* transform_data);

///@}
/**
 * DEPRECATED
 * @deprecated NOT USED
 */
void fclaw3d_patch_create_user_data(struct fclaw_global* glob,
                                    struct fclaw_patch* patch);

/**
 * DEPRECATED
 * @deprecated NOT USED
 */
void fclaw3d_patch_destroy_user_data(struct fclaw_global* glob,
                                     struct fclaw_patch* patch);

/* ------------------------------------------------------------------------------------ */
///                         @name Transform Functions
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Initialize the transform data for a patch
 * 
 * @param[in] glob the global context
 * @param[in] patch the patch context
 * @param[in] blockno the block number, -1 if ghost patch
 * @param[in] patchno the patch number
 * @param[in,out] tdata the stransform data structure
 */
void fclaw3d_patch_transform_init_data(struct fclaw_global* glob,
                                       struct fclaw_patch* patch,
                                       int blockno, int patchno,
                                       struct fclaw_patch_transform_data *tdata);

/**
 * @brief Get the transform on a block face
 * 
 * @param[in] glob the global context
 * @param[in] faceno 
 * @param[in] rfaceno 
 * @param[out]  ftransform  This array holds 9 integers.
 *              [0,2]       The coordinate axis sequence of the origin face,
 *                          the first referring to the tangential and the second
 *                          to the normal.  A permutation of (0, 1).
 *              [3,5]       The coordinate axis sequence of the target face.
 *              [6,8]       Edge reversal flag for tangential axis (boolean);
 *                          face code in [0, 3] for the normal coordinate q:
 *                          0: q' = -q
 *                          1: q' = q + 1
 *                          2: q' = q - 1
 *                          3: q' = 2 - q
 *                          [8] & 4: Both patches are in the same block,
 *                                   the \a ftransform contents are ignored.
 */
void fclaw3d_patch_transform_blockface(struct fclaw_global* glob,
                                       int faceno, int rfaceno,
                                       int ftransform[]);

/**
 * @brief Get the transform for within a block (the identity transform)
 * 
 * @param[in] glob the global context
 * @param[out]  ftransform  This array holds 9 integers.
 *              [0,2]       The coordinate axis sequence of the origin face,
 *                          the first referring to the tangential and the second
 *                          to the normal.  A permutation of (0, 1).
 *              [3,5]       The coordinate axis sequence of the target face.
 *              [6,8]       Edge reversal flag for tangential axis (boolean);
 *                          face code in [0, 3] for the normal coordinate q:
 *                          0: q' = -q
 *                          1: q' = q + 1
 *                          2: q' = q - 1
 *                          3: q' = 2 - q
 *                          [8] & 4: Both patches are in the same block,
 *                                   the \a ftransform contents are ignored.
 */
void fclaw3d_patch_transform_blockface_intra(struct fclaw_global* glob, 
                                             int ftransform[]);
  
///@}
/* ------------------------------------------------------------------------------------ */
///                         @name Regridding Functions
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Tag a patch for refinement
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] this_blockno the block number
 * @param[in] this_patchno the patch number
 * @param[in] initflag true if in init phase
 * @return true if patch should be refined
 */
int fclaw3d_patch_tag4refinement(struct fclaw_global *glob,
                                 struct fclaw_patch *this_patch,
                                 int blockno, int patchno,
                                 int initflag);

/**
 * @brief Tag a patch for coarsening
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] this_blockno the block number
 * @param[in] this_patchno the patch number
 * @param[in] initflag true if in init phase
 * @return true if patch should be coarsened
 */
int fclaw3d_patch_tag4coarsening(struct fclaw_global *glob,
                                 struct fclaw_patch *fine_patches,
                                 int blockno,
                                 int patchno,
                                 int initflag);

/**
 * @brief Interpolates a set of patches from a coarse patch
 * 
 * @param[in] glob the global context
 * @param[in] coarse_patch the coarse patch context
 * @param[in,out] fine_patches the fine patch contexts
 * @param[in] blockno the block number
 * @param[in] coarse_patchno the patch number of the coarse patch
 * @param[in] fine0_patchno the patch number of the first fine patch
 */
void fclaw3d_patch_interpolate2fine(struct fclaw_global *glob,
                                    struct fclaw_patch* coarse_patch,
                                    struct fclaw_patch* fine_patches,
                                    int this_blockno, int coarse_patchno,
                                    int fine0_patchno);

/**
 * @brief Averages from a set of fine patches to a coarse patch
 * 
 * @param[in] glob the global context
 * @param[in] fine_patches the fine patch contexts
 * @param[in,out] coarse_patch the coarse patch context
 * @param[in] blockno the block number
 * @param[in] fine_patchno the patch number of the first fine patch
 * @param[in] coarse_patchno the patch number of the coarse patch
 */
void fclaw3d_patch_average2coarse(struct fclaw_global *glob,
                                  struct fclaw_patch *fine_patches,
                                  struct fclaw_patch *coarse_patch,
                                  int blockno, int fine0_patchno,
                                  int coarse_patchno);

///@}
/* ------------------------------------------------------------------------------------ */
///                         @name Parallel Ghost Patches
/* ------------------------------------------------------------------------------------ */
///@{

/** 
 * @brief Get the buffer size needed to pack a single patch's ghost data 
 * 
 * @param[in] glob the global context
 * @return the buffer size (in bytes)
 */
size_t fclaw3d_patch_ghost_packsize(struct fclaw_global* glob);

/**
 * @brief Allocates a buffer for the patch ghost data
 * 
 * @param[in] glob the global context
 * @param[out] q pointer to the allocated buffer
 */
void fclaw3d_patch_local_ghost_alloc(struct fclaw_global* glob,
                                     void** q);

/**
 * @brief Frees a buffer for the patch ghost data
 * 
 * @param[in] glob the global context
 * @param[out] q pointer to the buffer to free
 */
void fclaw3d_patch_local_ghost_free(struct fclaw_global* glob,
                                    void **q);

/**
 * @brief Packs the patch ghost data into a buffer
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in,out] patch_data the buffer
 * @param[in] time_interp true if ghost filling for time interpolated level (non-global update)
 */
void fclaw3d_patch_local_ghost_pack(struct fclaw_global *glob,
                                    struct fclaw_patch *this_patch,
                                    void *patch_data,
                                    int time_interp);

/**
 * @brief Builds a new ghost patch
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in] build_mode the build mode
 */
void fclaw3d_patch_remote_ghost_build(struct fclaw_global *glob,
                                      struct fclaw_patch *this_patch,
                                      int blockno,
                                      int patchno,
                                      fclaw_build_mode_t build_mode);

/**
 * @brief Unpacks a ghost patch from a buffer
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number 
 * @param[in] patchno the patch number
 * @param[in] qdata the buffer to unpack from
 * @param[in] time_interp true if ghost filling for time interpolated level (non-global update)
 */
void fclaw3d_patch_remote_ghost_unpack(struct fclaw_global* glob,
                                       struct fclaw_patch* this_patch,
                                       int blockno, int patchno,
                                       void *qdata, int time_interp);


/**
 * @brief Frees memory used by a ghost patch
 * 
 * @param[in] glob the global context
 * @param[in,out] ghost_patch the patch context
 */
void fclaw3d_patch_remote_ghost_delete(struct fclaw_global *glob,
                                       struct fclaw_patch *ghost_patch);

///@}
/* ------------------------------------------------------------------------------------ */
///                          @name Parallel Partitioning
/* ------------------------------------------------------------------------------------ */
///@{


/**
 * @brief Packs a patch into a buffer
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[out] pack_data_here the buffer
 */
void fclaw3d_patch_partition_pack(struct fclaw_global *glob,
                                  struct fclaw_patch *this_patch,
                                  int blockno,
                                  int patchno,
                                  void *pack_data_here);

/**
 * @brief Unpacks a patch from a buffer
 * 
 * @param[in] glob the global context (contains the old domain)
 * @param[in] new_domain the new domain
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in] packed_data the buffer
 */
void fclaw3d_patch_partition_unpack(struct fclaw_global *glob,
                                    struct fclaw_domain *new_domain,  
                                    struct fclaw_patch *this_patch,
                                    int blockno,
                                    int patchno,
                                    void *packed_data);

/**
 * @brief Gets the buffer size (in bytes) needed to pack a patch
 * 
 * @param[in] glob the global context
 * @return size_t the size of buffer needed
 */
size_t fclaw3d_patch_partition_packsize(struct fclaw_global* glob);


///@}
/* ------------------------------------------------------------------------------------ */
///                         @name Time Syncing Functions
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Adds fine grid corrections to coarse grid.  
 * 
 * @param[in] glob the global context
 * @param[in,out] coarse_patch the coarse patch context
 * @param[in] fine_patch the fine patch context
 * @param[in] coarse_blockno the block number of the coarse patch
 * @param[in] fine_blockno the block number of the fine patch
 * @param[in] coarse_patchno the patch number of the coarse patch
 * @param[in] idir the direction of the interface 0 for bottom/top 
 *            1 for left/right
 * @param[in] igrid the index of the fine grid in the child array
 * @param[in] iface_coarse the interface on the coarse patch
 * @param[in] time_interp true if ghost filling for time interpolated level (non-global update)
 * @param[in] transform_data the transform for the neighbor's coordinates
 */
void fclaw3d_patch_time_sync_f2c(struct fclaw_global* glob,
                                 struct fclaw_patch *coarse_patch,
                                 struct fclaw_patch *fine_patch,
                                 int coarse_blockno, int fine_blockno,
                                 int coarse_patchno, 
                                 int idir,
                                 int igrid,
                                 int iface_coarse,
                                 int time_interp,
                                 struct fclaw_patch_transform_data* transform_data);

/**
 * @brief Adds corrections to patches that are at the same levle and are at block boundaries.
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch this patch
 * @param[in] neighbor_patch the neighbor patch
 * @param[in] this_iface the interface that the neighbor patch is on
 * @param[in] idir the direction of the interface 0 for bottom/top 
 *            1 for left/right
 * @param[in] transform_data the transform for the neighbor's coordinates
 */
void fclaw3d_patch_time_sync_samesize(struct fclaw_global* glob,
                                      struct fclaw_patch *this_patch,
                                      struct fclaw_patch *neighbor_patch,
                                      int this_iface, int idir,
                                      struct fclaw_patch_transform_data *transform_data);

/**
 * @brief Resets conservation data
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] coarse_level the the level of the coarse patch
 * @param[in] reset_mode the reset mode ::fclaw3d_time_sync_type
 */
void fclaw3d_patch_time_sync_reset(struct fclaw_global* glob,
                                   struct fclaw_patch* this_patch,
                                   int coarse_level,
                                   int reset_mode);


///@}
/* ------------------------------------------------------------------------------------ */
///                            @name Misc Access Functions
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Get the global_num, local_num, and level of a patch
 * 
 * @param[in] domain 
 * @param[in] patch 
 * @param[in] blockno 
 * @param[in] patchno 
 * @param[out] global_num the global patch number
 * @param[out] local_num the local patch number
 * @param[out] level the level that the patch is on
 */
void fclaw3d_patch_get_info(struct fclaw_domain * domain,
                            struct fclaw_patch * patch,
                            int blockno, int patchno,
                            int *global_num, int* local_num, 
                            int *level);

/**
 * @brief Get the block number, patch number, global_num, and level of a patch
 * 
 * @param[in] domain the domain
 * @param[in] this_patch the patch context
 * @param[out] blockno the block number
 * @param[out] patchno the patch number
 * @param[out] global_num the global patch number
 * @param[out] level the level
 */
/* I don't completely trust this routine */
void fclaw3d_patch_get_info2(struct fclaw_domain * domain,
                             struct fclaw_patch * this_patch,
                             int *blockno, int *patchno,
                             int *global_num, int *level);
/**
 * @brief Get the user patch pointer
 * 
 * @param patch the patch context
 * @return void* the pointer
 */
void* fclaw3d_patch_get_user_patch(struct fclaw_patch* patch);

/**
 * @brief Get the patch data
 * 
 * @param patch the patch context
 * @return struct fclaw3d_patch_data* pointer to the patch data
 */
struct fclaw_patch_data* fclaw3d_patch_get_patch_data(struct fclaw_patch* patch);

/**
 * @brief Get the user data pointer
 * 
 * @param glob the global context
 * @param this_patch the patch context
 * @return void* the user data pointer
 */
void* fclaw3d_patch_get_user_data(struct fclaw_global* glob,
                                  struct fclaw_patch* this_patch);


/**
 * @brief Get the metric patch
 * 
 * @param glob the global context
 * @param patch the patch context
 * @return void* pointer to the metric patch
 */
void* fclaw3d_patch_metric_patch(struct fclaw_global* glob,
                                 struct fclaw_patch *patch);

/**
 * @brief Get the block number
 * 
 * @param this_patch the patch context
 * @return int the block number
 */
int fclaw3d_patch_get_blockno(struct fclaw_patch* this_patch);

/**
 * @brief Get the patch number
 * 
 * @param this_patch the patch context
 * @return int the patch number
 */
int fclaw3d_patch_get_patchno(struct fclaw_patch* this_patch);

///@}
/* ------------------------------------------------------------------------------------ */
///                                 @name Misc User Data
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Get the user data pointer of a patch
 * 
 * @param glob the global context
 * @param this_patch the patch context
 * @return void* the user data pointer
 */
void* fclaw3d_patch_user_data(struct fclaw_global* glob,
                              struct fclaw_patch* this_patch);

/**
 * @brief Set the user data pointer of a patch
 * 
 * @param glob the global context
 * @param this_patch the patch context
 * @param user the user data pointer
 */
void fclaw3d_patch_set_user_data(struct fclaw_global* glob,
                                 struct fclaw_patch* this_patch, 
                                 void* user);


///@}
/* ------------------------------------------------------------------------------------ */
///                       @name Misc Functions (mostly internal)
/* ------------------------------------------------------------------------------------ */
///@{

/**
 * @brief Returns true if patch lies on a parallel boundary
 * 
 * @param patch the patch context
 * @return int true if on parallel boundary
 */
int fclaw3d_patch_on_parallel_boundary (const struct fclaw_patch * patch);


/**
 * @brief Set the face type for a patch
 * 
 * @param patch the patch context
 * @param iface the interface
 * @param face_type the face type
 */
void fclaw3d_patch_set_face_type(struct fclaw_patch *patch, int iface,
                                 fclaw_patch_relation_t face_type);

/**
 * @brief Set the edge type for a patch
 * 
 * @param patch the patch context
 * @param iedge the edge
 * @param edge_type the edge type
 */
void fclaw3d_patch_set_edge_type(fclaw_patch_t *patch,int iedge,
								   fclaw_patch_relation_t edge_type);

/**
 * @brief Set the corner type for a patch
 * 
 * @param patch the patch context
 * @param icorner the corner
 * @param corner_type the corner type
 */
void fclaw3d_patch_set_corner_type(struct fclaw_patch *patch, int icorner,
                                   fclaw_patch_relation_t corner_type);

/**
 * @brief Set the missing corner of a patch
 * 
 * @param patch the patch context
 * @param icorner the missing corner
 */
void fclaw3d_patch_set_missing_corner(struct fclaw_patch *patch, int icorner);

/**
 * @brief Get the face type of a patch
 * 
 * @param patch the patch context
 * @param iface the face
 * @return fclaw3d_patch_relation_t the face type
 */
fclaw_patch_relation_t fclaw3d_patch_get_face_type(struct fclaw_patch* patch,
                                                        int iface);

/**
 * @brief Get the edge type of a patch
 * 
 * @param patch the patch context
 * @param iface the edge
 * @return fclaw3d_patch_relation_t the edge type
 */
fclaw_patch_relation_t fclaw2d_patch_get_edge_type(fclaw_patch_t* patch,
													   int iedge);

/**
 * @brief Get the corner type of a patch
 * 
 * @param patch the patch context
 * @param icorner the corner
 * @return fclaw3d_patch_relation_t the patch relation
 */
fclaw_patch_relation_t fclaw3d_patch_get_corner_type(struct fclaw_patch* patch,
                                                          int icorner);

/**
 * @brief Returns true if a corner is missing
 * 
 * @param patch the patch context
 * @param icorner the corner
 * @return int true if a corner is missing
 */
int fclaw3d_patch_corner_is_missing(struct fclaw_patch* patch,
                                    int icorner);

/**
 * @brief Set the neighbor relation data for a patch 
 * 
 * @param patch the patch context
 */
void fclaw3d_patch_neighbors_set(fclaw_patch_t* patch);

/**
 * @brief Reset the neighbor relation data for a patch
 * 
 * @param patch the patch context
 */
void fclaw3d_patch_neighbors_reset(struct fclaw_patch* patch);

/**
 * @brief Returns true if the patch neighbor information is set
 * 
 * @param patch the patch context
 * @return int true if the patch neighbor information is set
 */
int fclaw3d_patch_neighbor_type_set(struct fclaw_patch* patch);

/**
 * @brief Returns true if a patch has finer neighbors
 * 
 * @param patch the patch context
 * @return int true if the patch has finer neighbors
 */
int fclaw3d_patch_has_finegrid_neighbors(struct fclaw_patch *patch);

/**
 * @brief Returns true if the patch is on a coarse fine interface
 * 
 * @param patch the patch context
 * @return int true if the patch is on a coarse fine interface
 */
int fclaw3d_patch_on_coarsefine_interface(struct fclaw_patch *patch);

/**
 * @brief Get the block corner count array for a patch
 * 
 * @param glob the global context
 * @param this_patch the patch context
 * @return int* the array with the block corner count (the number of blocks that meet) for each corner
 */
int* fclaw3d_patch_block_corner_count(struct fclaw_global *glob,
                                      struct fclaw_patch* this_patch);

/**
 * @brief Set the block corner count for a corner
 * 
 * @param glob the global context
 * @param this_patch the patch context
 * @param icorner the corner to set
 * @param block_corner_count the block corner count (the number of blocks that meet)
 */
void fclaw3d_patch_set_block_corner_count(struct fclaw_global *glob,
                                          struct fclaw_patch* this_patch,
                                          int icorner, int block_corner_count);

///@}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
