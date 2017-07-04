/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW2D_PATCH_H
#define FCLAW2D_PATCH_H

#include <forestclaw2d.h>  /* Contains definition of patch-iterator callback */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fclaw2d_patch_vtable  fclaw2d_patch_vtable_t;
typedef struct fclaw2d_patch_data    fclaw2d_patch_data_t;

typedef enum
{
    FCLAW2D_BUILD_FOR_GHOST_AREA_COMPUTED = 0,
    FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED,
    FCLAW2D_BUILD_FOR_UPDATE,
    FCLAW2D_BUILD_CUSTOM
} fclaw2d_build_mode_t;


/* The user patch (clawpatch, for example) is stored as a 'user_patch', below. */
struct fclaw2d_patch_data
{
    fclaw2d_patch_relation_t face_neighbors[4];
    fclaw2d_patch_relation_t corner_neighbors[4];
    int corners[4];
    int block_corner_count[4];
    int on_coarsefine_interface;
    int has_finegrid_neighbors;
    int neighbors_set;

    void *user_patch; /* Start of attempt to "virtualize" the user patch. */
};

struct fclaw2d_global;
struct fclaw2d_domain;
struct fclaw2d_patch;            
struct fclaw2d_transform_data;


/* ------------------------------ Creating/deleting patches --------------------------- */

#if 0
void fclaw2d_patch_data_new(struct fclaw2d_global* glob,
                            struct fclaw2d_patch* this_patch);
#endif                            

void fclaw2d_patch_data_delete(struct fclaw2d_global *glob,
                               struct fclaw2d_patch *patch);

void fclaw2d_patch_build(struct fclaw2d_global *glob,
                         struct fclaw2d_patch *this_patch,
                         int blockno,
                         int patchno,
                         void *user);

void fclaw2d_patch_build_from_fine(struct fclaw2d_global *glob,
                                   struct fclaw2d_patch *fine_patches,
                                   struct fclaw2d_patch *coarse_patch,
                                   int blockno,
                                   int coarse_patchno,
                                   int fine0_patchno,
                                   fclaw2d_build_mode_t build_mode);

/* ---------------------------- Solver specific functions ----------------------------- */

void fclaw2d_patch_initialize(struct fclaw2d_global *glob,
                              struct fclaw2d_patch *this_patch,
                              int this_block_idx,
                              int this_patch_idx);

void fclaw2d_patch_physical_bc(struct fclaw2d_global *glob,
                               struct fclaw2d_patch *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               double t,
                               double dt,
                               int *intersects_bc,
                               int time_interp);

double fclaw2d_patch_single_step_update(struct fclaw2d_global *glob,
                                        struct fclaw2d_patch *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        double t,
                                        double dt);

/* -------------------------------- time stepping ------------------------------------- */

void fclaw2d_patch_restore_step(struct fclaw2d_global* glob,
                                struct fclaw2d_patch* this_patch);

void fclaw2d_patch_save_step(struct fclaw2d_global* glob,
                             struct fclaw2d_patch* this_patch);


void fclaw2d_patch_setup_timeinterp(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch *this_patch,
                                    double alpha);


/* ---------------------------- Ghost filling - patch specific ------------------------ */

void fclaw2d_patch_copy_face(struct fclaw2d_global* glob,
                             struct fclaw2d_patch *this_patch,
                             struct fclaw2d_patch *neighbor_patch,
                             int iface,
                             int time_interp,
                             struct fclaw2d_transform_data *transform_data);

void fclaw2d_patch_average_face(struct fclaw2d_global* glob,
                                struct fclaw2d_patch *coarse_patch,
                                struct fclaw2d_patch *fine_patch,
                                int idir,
                                int iface_coarse,
                                int RefineFactor,
                                int refratio,
                                int time_interp,
                                int igrid,
                                struct fclaw2d_transform_data* transform_data);

void fclaw2d_patch_interpolate_face(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch *coarse_patch,
                                    struct fclaw2d_patch *fine_patch,
                                    int idir,
                                    int iside,
                                    int RefineFactor,
                                    int refratio,
                                    int time_interp,
                                    int igrid,
                                    struct fclaw2d_transform_data* transform_data);


void fclaw2d_patch_copy_corner(struct fclaw2d_global* glob,
                               struct fclaw2d_patch *this_patch,
                               struct fclaw2d_patch *corner_patch,
                               int icorner,
                               int time_interp,
                               struct fclaw2d_transform_data *transform_data);

void fclaw2d_patch_average_corner(struct fclaw2d_global* glob,
                                  struct fclaw2d_patch *coarse_patch,
                                  struct fclaw2d_patch *fine_patch,
                                  int coarse_corner,
                                  int refratio,
                                  int time_interp,
                                  struct fclaw2d_transform_data* transform_data);



void fclaw2d_patch_interpolate_corner(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch* coarse_patch,
                                      struct fclaw2d_patch* fine_patch,
                                      int coarse_corner,
                                      int refratio,
                                      int time_interp,
                                      struct fclaw2d_transform_data* transform_data);


/* ------------------------------- Regridding functions ------------------------------- */

int fclaw2d_patch_tag4refinement(struct fclaw2d_global *glob,
                                 struct fclaw2d_patch *this_patch,
                                 int blockno, int patchno,
                                 int initflag);

int fclaw2d_patch_tag4coarsening(struct fclaw2d_global *glob,
                                 struct fclaw2d_patch *fine_patches,
                                 int blockno,
                                 int patchno);

void fclaw2d_patch_interpolate2fine(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch* coarse_patch,
                                    struct fclaw2d_patch* fine_patches,
                                    int this_blockno, int coarse_patchno,
                                    int fine0_patchno);

void fclaw2d_patch_average2coarse(struct fclaw2d_global *glob,
                                  struct fclaw2d_patch *fine_patches,
                                  struct fclaw2d_patch *coarse_patch,
                                  int blockno, int fine0_patchno,
                                  int coarse_patchno);

/* ----------------------------- Parallel ghost patches ------------------------------- */

size_t fclaw2d_patch_ghost_packsize(struct fclaw2d_global* glob);


void fclaw2d_patch_local_ghost_alloc(struct fclaw2d_global* glob,
                                     void** q);

void fclaw2d_patch_local_ghost_free(struct fclaw2d_global* glob,
                                    void **q);

void fclaw2d_patch_local_ghost_pack(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch *this_patch,
                                    void *patch_data,
                                    int time_interp);

void fclaw2d_patch_remote_ghost_build(struct fclaw2d_global *glob,
                                      struct fclaw2d_patch *this_patch,
                                      int blockno,
                                      int patchno,
                                      void *user);

void fclaw2d_patch_remote_ghost_unpack(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch* this_patch,
                                       int this_block_idx, int this_patch_idx,
                                       void *qdata, int time_interp);


void fclaw2d_patch_remote_ghost_delete(struct fclaw2d_global *glob,
                                       struct fclaw2d_patch *ghost_patch);

/* -------------------------------- Parallel partitioning ----------------------------- */


void fclaw2d_patch_partition_pack(struct fclaw2d_global *glob,
                                  struct fclaw2d_patch *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx,
                                  void *pack_data_here);

void fclaw2d_patch_partition_unpack(struct fclaw2d_global *glob,  /* contains old domain */
                                    struct fclaw2d_domain *new_domain,  
                                    struct fclaw2d_patch *this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    void *packed_data);

size_t fclaw2d_patch_partition_packsize(struct fclaw2d_global* glob);


/* ------------------------------ Misc access functions ------------------------------- */
void fclaw2d_patch_get_info(struct fclaw2d_domain * domain,
                            struct fclaw2d_patch * this_patch,
                            int this_block_idx, int this_patch_idx,
                            int *global_num, int *level);


void*
fclaw2d_patch_get_user_patch(fclaw2d_patch_t* patch);

struct fclaw2d_patch_data*
fclaw2d_patch_get_user_data(struct fclaw2d_patch* patch);


void* fclaw2d_patch_get_user_patch(struct fclaw2d_patch* patch);


/* ---------------------- Creating/deleting patches (typedefs) ------------------------ */

typedef void* (*fclaw2d_patch_new_t)();

typedef void (*fclaw2d_patch_delete_t)(void *user_patch);

typedef void (*fclaw2d_patch_build_t)(struct fclaw2d_global *glob,
                                      struct fclaw2d_patch *this_patch,
                                      int blockno,
                                      int patchno,
                                      void *user);

typedef void (*fclaw2d_patch_build_from_fine_t)(struct fclaw2d_global *glob,
                                                struct fclaw2d_patch *fine_patches,
                                                struct fclaw2d_patch *coarse_patch,
                                                int blockno,
                                                int coarse_patchno,
                                                int fine0_patchno,
                                                fclaw2d_build_mode_t build_mode);

typedef void (*fclaw2d_patch_setup_t)(struct fclaw2d_global *glob,
                                      struct fclaw2d_patch *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx);

/* --------------------- Solver specific functions (typedefs) ------------------------- */

typedef void (*fclaw2d_patch_initialize_t)(struct fclaw2d_global *glob,
                                           struct fclaw2d_patch *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx);

typedef void (*fclaw2d_patch_physical_bc_t)(struct fclaw2d_global *glob,
                                            struct fclaw2d_patch *this_patch,
                                            int this_block_idx,
                                            int this_patch_idx,
                                            double t,
                                            double dt,
                                            int *intersects_bc,
                                            int time_interp);

typedef double (*fclaw2d_patch_single_step_update_t)(struct fclaw2d_global *glob,
                                                     struct fclaw2d_patch *this_patch,
                                                     int this_block_idx,
                                                     int this_patch_idx,
                                                     double t,
                                                     double dt);

/* ----------------------------- Time stepping (typedefs) ----------------------------- */

typedef void (*fclaw2d_patch_setup_timeinterp_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *this_patch,
                                                 double alpha);

typedef void (*fclaw2d_patch_restore_step_t)(struct fclaw2d_global *glob,
                                             struct fclaw2d_patch* this_patch);

typedef void (*fclaw2d_patch_save_step_t)(struct fclaw2d_global *glob,
                                          struct fclaw2d_patch* this_patch);


/* --------------------- Ghost filling - patch specific (typedefs) -------------------- */

typedef void (*fclaw2d_patch_copy_face_t)(struct fclaw2d_global* glob,
                                          struct fclaw2d_patch *this_patch,
                                          struct fclaw2d_patch *neighbor_patch,
                                          int iface,
                                          int time_interp,
                                          struct fclaw2d_transform_data *transform_data);

typedef void (*fclaw2d_patch_average_face_t)(struct fclaw2d_global* glob,
                                             struct fclaw2d_patch *coarse_patch,
                                             struct fclaw2d_patch *fine_patch,
                                             int idir,
                                             int iface_coarse,
                                             int RefineFactor,
                                             int refratio,
                                             int time_interp,
                                             int igrid,
                                             struct fclaw2d_transform_data* transform_data);

typedef void (*fclaw2d_patch_interpolate_face_t)(struct fclaw2d_global* glob,
                                                 struct fclaw2d_patch *coarse_patch,
                                                 struct fclaw2d_patch *fine_patch,
                                                 int idir,
                                                 int iside,
                                                 int RefineFactor,
                                                 int refratio,
                                                 int a_time_interp,
                                                 int igrid,
                                                 struct fclaw2d_transform_data* transform_data);

typedef void (*fclaw2d_patch_copy_corner_t)(struct fclaw2d_global* glob,
                                            struct fclaw2d_patch *this_patch,
                                            struct fclaw2d_patch *corner_patch,
                                            int icorner,
                                            int time_interp,
                                            struct fclaw2d_transform_data *transform_data);

typedef void (*fclaw2d_patch_average_corner_t)(struct fclaw2d_global* glob,
                                               struct fclaw2d_patch *coarse_patch,
                                               struct fclaw2d_patch *fine_patch,
                                               int coarse_corner,
                                               int refratio,
                                               int time_interp,
                                               struct fclaw2d_transform_data* transform_data);

typedef void (*fclaw2d_patch_interpolate_corner_t)(struct fclaw2d_global* glob,
                                                   struct fclaw2d_patch* coarse_patch,
                                                   struct fclaw2d_patch* fine_patch,
                                                   int coarse_corner,
                                                   int refratio,
                                                   int a_time_interp,
                                                   struct fclaw2d_transform_data* transform_data);

/* ------------------------- Regridding functions (typedefs) -------------------------- */

typedef int (*fclaw2d_patch_tag4refinement_t)(struct fclaw2d_global *glob,
                                              struct fclaw2d_patch *this_patch,
                                              int this_block_idx, int this_patch_idx,
                                              int initflag);

typedef int (*fclaw2d_patch_tag4coarsening_t)(struct fclaw2d_global *glob,
                                               struct fclaw2d_patch *this_patch,
                                               int this_blockno,
                                               int this_patchno);

typedef void (*fclaw2d_patch_interpolate2fine_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *coarse_patch,
                                                 struct fclaw2d_patch* fine_patches,
                                                 int this_blockno, int coarse_patchno,
                                                 int fine_patchno);

typedef void (*fclaw2d_patch_average2coarse_t)(struct fclaw2d_global *glob,
                                               struct fclaw2d_patch *fine_siblings,
                                               struct fclaw2d_patch *coarse_patch,
                                               int blockno, int fine_patchno,
                                               int coarse_patchno);

/* -------------------------- Parallel ghost patches (typedefs) ----------------------- */

typedef size_t (*fclaw2d_patch_ghost_packsize_t)(struct fclaw2d_global* glob);

typedef void (*fclaw2d_patch_local_ghost_pack_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *this_patch,
                                                 void *patch_data,
                                                 int time_interp);

typedef void (*fclaw2d_patch_local_ghost_alloc_t)(struct fclaw2d_global* glob,
                                                 void** q);

typedef void (*fclaw2d_patch_local_ghost_free_t)(struct fclaw2d_global* glob,
                                                 void **q);

typedef void (*fclaw2d_patch_remote_ghost_build_t)(struct fclaw2d_global *glob,
                                                   struct fclaw2d_patch *this_patch,
                                                   int blockno,
                                                   int patchno,
                                                   void *user);

typedef void (*fclaw2d_patch_remote_ghost_setup_t)(struct fclaw2d_global *glob,
                                                   struct fclaw2d_patch *this_patch,
                                                   int this_block_idx,
                                                   int this_patch_idx);


typedef void (*fclaw2d_patch_remote_ghost_unpack_t)(struct fclaw2d_global *glob,
                                                    struct fclaw2d_patch* this_patch,
                                                    int this_block_idx, int this_patch_idx,
                                                    void *qdata, int time_interp);

typedef void (*fclaw2d_patch_remote_ghost_delete_t)(void *user_patch);



/* ----------------------------- Partitioning (typedefs) ------------------------------ */

/* Returns size, in bytes, i.e. psize = mx*my*sizeof(double) */ 
typedef size_t (*fclaw2d_patch_partition_packsize_t)(struct fclaw2d_global* glob);


typedef void (*fclaw2d_patch_partition_pack_t)(struct fclaw2d_global *glob,
                                               struct fclaw2d_patch *this_patch,
                                               int this_block_idx,
                                               int this_patch_idx,
                                               void *pack_data_here);

typedef void (*fclaw2d_patch_partition_unpack_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_domain *new_domain,
                                                 struct fclaw2d_patch *this_patch,
                                                 int this_block_idx,
                                                 int this_patch_idx,
                                                 void *unpack_data_from_here);

/* ----------------------------------- Virtual table  --------------------------------- */
struct fclaw2d_patch_vtable
{

    /* Solver functions */
    fclaw2d_patch_initialize_t            initialize;
    fclaw2d_patch_physical_bc_t           physical_bc;
    fclaw2d_patch_single_step_update_t    single_step_update;

    /* Creating/deleting/building patches */
    fclaw2d_patch_new_t                   patch_new;
    fclaw2d_patch_delete_t                patch_delete;
    fclaw2d_patch_build_t                 build;
    fclaw2d_patch_build_from_fine_t       build_from_fine;
    fclaw2d_patch_setup_t                 setup;

    /* Time stepping */
    fclaw2d_patch_restore_step_t          restore_step;
    fclaw2d_patch_save_step_t             save_step;
    fclaw2d_patch_setup_timeinterp_t      setup_timeinterp;


    /* regridding functions */
    fclaw2d_patch_tag4refinement_t        tag4refinement;
    fclaw2d_patch_tag4coarsening_t        tag4coarsening;
    fclaw2d_patch_average2coarse_t        average2coarse;
    fclaw2d_patch_interpolate2fine_t      interpolate2fine;

    /* ghost filling functions */
    fclaw2d_patch_copy_face_t             copy_face;
    fclaw2d_patch_average_face_t          average_face;
    fclaw2d_patch_interpolate_face_t      interpolate_face;

    fclaw2d_patch_copy_corner_t           copy_corner;
    fclaw2d_patch_average_corner_t        average_corner;
    fclaw2d_patch_interpolate_corner_t    interpolate_corner;

    /* Ghost packing functions (for parallel use) */
    fclaw2d_patch_ghost_packsize_t        ghost_packsize;
    fclaw2d_patch_local_ghost_pack_t      local_ghost_pack;
    fclaw2d_patch_local_ghost_alloc_t     local_ghost_alloc;
    fclaw2d_patch_local_ghost_free_t      local_ghost_free;

    fclaw2d_patch_remote_ghost_build_t    remote_ghost_build;
    fclaw2d_patch_remote_ghost_setup_t    remote_ghost_setup;   /* Remote ghost patches */
    fclaw2d_patch_remote_ghost_unpack_t   remote_ghost_unpack;
    fclaw2d_patch_remote_ghost_delete_t   remote_ghost_delete;  /* Delete remote ghosts */

    /* Parallel load balancing (partitioning) */
    fclaw2d_patch_partition_pack_t         partition_pack;
    fclaw2d_patch_partition_unpack_t       partition_unpack;
    fclaw2d_patch_partition_packsize_t     partition_packsize;

    int is_set;
};

fclaw2d_patch_vtable_t* fclaw2d_patch_vt();

void fclaw2d_patch_vtable_initialize();


/* ----------------------------  Other access functions ------------------------------- */

typedef void (*fclaw2d_patch_iterator_t) (struct fclaw2d_global * glob, int level,
                                          fclaw2d_patch_callback_t pcb, void *user);

int fclaw2d_patch_on_parallel_boundary (const struct fclaw2d_patch * patch);


void fclaw2d_patch_set_face_type(struct fclaw2d_patch *patch, int iface,
                                 fclaw2d_patch_relation_t face_type);

void fclaw2d_patch_set_corner_type(struct fclaw2d_patch *patch, int icorner,
                                   fclaw2d_patch_relation_t corner_type);

void fclaw2d_patch_set_missing_corner(struct fclaw2d_patch *patch, int icorner);

fclaw2d_patch_relation_t fclaw2d_patch_get_face_type(struct fclaw2d_patch* patch,
                                                        int iface);
fclaw2d_patch_relation_t fclaw2d_patch_get_corner_type(struct fclaw2d_patch* patch,
                                                          int icorner);

int fclaw2d_patch_corner_is_missing(struct fclaw2d_patch* patch,
                                    int icorner);

void fclaw2d_patch_neighbors_set(struct fclaw2d_patch* patch);

void fclaw2d_patch_neighbors_reset(struct fclaw2d_patch* patch);

int fclaw2d_patch_has_finegrid_neighbors(struct fclaw2d_patch *patch);

int fclaw2d_patch_on_coarsefine_interface(struct fclaw2d_patch *patch);

int* fclaw2d_patch_block_corner_count(struct fclaw2d_global *glob,
                                      struct fclaw2d_patch* this_patch);

void fclaw2d_patch_set_block_corner_count(struct fclaw2d_global *glob,
                                          struct fclaw2d_patch* this_patch,
                                          int icorner, int block_corner_count);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
