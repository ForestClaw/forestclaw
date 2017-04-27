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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_transform.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef enum
{
    FCLAW2D_BUILD_FOR_GHOST_AREA_COMPUTED = 0,
    FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED,
    FCLAW2D_BUILD_FOR_UPDATE,
    FCLAW2D_BUILD_COSTOM
} fclaw2d_build_mode_t;

typedef void (*fclaw2d_patch_iterator_t) (fclaw2d_global_t * glob, int level,
                                          fclaw2d_patch_callback_t pcb, void *user);


/* Opaque pointer */
typedef struct fclaw2d_patch_data fclaw2d_patch_data_t;

void fclaw2d_patch_data_new(fclaw2d_global_t* glob,
                            fclaw2d_patch_t* this_patch);

void fclaw2d_patch_data_delete(fclaw2d_global_t *glob,
                               fclaw2d_patch_t *patch);

void fclaw2d_patch_delete_remote_ghost(fclaw2d_global_t *glob,
                                       fclaw2d_patch_t *ghost_patch);

void fclaw2d_patch_initialize(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx);

void fclaw2d_patch_write_header(fclaw2d_global_t* glob,
                                int iframe);


void fclaw2d_patch_physical_bc(fclaw2d_global_t *glob,
                               fclaw2d_patch_t *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               double t,
                               double dt,
                               fclaw_bool *intersects_bc,
                               fclaw_bool time_interp);


struct fclaw2d_patch_data*
fclaw2d_patch_get_data(fclaw2d_patch_t* patch);


void fclaw2d_domain_iterate_level_mthread (fclaw2d_domain_t * domain, int level,
                                           fclaw2d_patch_callback_t pcb, void *user);

void* fclaw2d_patch_get_user_patch(fclaw2d_patch_t* patch);

/* --------------------------------------------------------------
   Routines that return information about connectivity.
   This information is obtained with each new regridding.
   ------------------------------------------------------------ */

int fclaw2d_patch_on_parallel_boundary (const fclaw2d_patch_t * patch);


void fclaw2d_patch_set_face_type(fclaw2d_patch_t *patch,int iface,
                                 fclaw2d_patch_relation_t face_type);

void fclaw2d_patch_set_corner_type(fclaw2d_patch_t *patch,int icorner,
                                   fclaw2d_patch_relation_t corner_type);

void fclaw2d_patch_set_missing_corner(fclaw2d_patch_t *patch,int icorner);

fclaw2d_patch_relation_t fclaw2d_patch_get_face_type(fclaw2d_patch_t* patch,
                                                        int iface);
fclaw2d_patch_relation_t fclaw2d_patch_get_corner_type(fclaw2d_patch_t* patch,
                                                          int icorner);

int fclaw2d_patch_corner_is_missing(fclaw2d_patch_t* patch,
                                    int icorner);

void fclaw2d_patch_neighbors_set(fclaw2d_patch_t* patch);

void fclaw2d_patch_neighbors_reset(fclaw2d_patch_t* patch);

int fclaw2d_patch_has_finegrid_neighbors(fclaw2d_patch_t *patch);

int fclaw2d_patch_on_coarsefine_interface(fclaw2d_patch_t *patch);

int* fclaw2d_patch_block_corner_count(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t* this_patch);

void fclaw2d_patch_set_block_corner_count(fclaw2d_global_t *glob,
                                          fclaw2d_patch_t* this_patch,
                                          int icorner, int block_corner_count);

void fclaw2d_patch_pack_local_ghost(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    double *patch_data,
                                    int time_interp);

void fclaw2d_patch_unpack_remote_ghost(fclaw2d_global_t* glob,
                                       fclaw2d_patch_t* this_patch,
                                       int this_block_idx, int this_patch_idx,
                                       double *qdata, fclaw_bool time_interp);

void fclaw2d_patch_build_remote_ghost(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno,
                                      void *user);

size_t fclaw2d_patch_ghost_packsize(fclaw2d_global_t* glob);

void fclaw2d_patch_alloc_local_ghost(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch,
                                     void** q);

void fclaw2d_patch_free_local_ghost(fclaw2d_global_t* glob,
                                    void **q);


void cb_fclaw2d_patch_partition_pack(fclaw2d_domain_t *domain,
                                     fclaw2d_patch_t *this_patch,
                                     int this_block_idx,
                                     int this_patch_idx,
                                     void *user);

void fclaw2d_patch_partition_unpack(fclaw2d_global_t *glob,
                                    fclaw2d_domain_t *new_domain,
                                    fclaw2d_patch_t *this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    void *user);

size_t fclaw2d_patch_partition_packsize(fclaw2d_global_t* glob);

void fclaw2d_patch_build(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int blockno,
                         int patchno,
                         void *user);

void fclaw2d_patch_build_from_fine(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *fine_patches,
                                   fclaw2d_patch_t *coarse_patch,
                                   int blockno,
                                   int coarse_patchno,
                                   int fine0_patchno,
                                   fclaw2d_build_mode_t build_mode);

void fclaw2d_patch_restore_step(fclaw2d_global_t* glob,
                                fclaw2d_patch_t* this_patch);

void fclaw2d_patch_save_step(fclaw2d_global_t* glob,
                             fclaw2d_patch_t* this_patch);

void fclaw2d_patch_interpolate_face(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t *coarse_patch,
                                    fclaw2d_patch_t *fine_patch,
                                    int idir,
                                    int iside,
                                    int p4est_refineFactor,
                                    int refratio,
                                    fclaw_bool time_interp,
                                    int igrid,
                                    fclaw2d_transform_data_t* transform_data);

void fclaw2d_patch_average_face(fclaw2d_global_t* glob,
                                fclaw2d_patch_t *coarse_patch,
                                fclaw2d_patch_t *fine_patch,
                                int idir,
                                int iface_coarse,
                                int p4est_refineFactor,
                                int refratio,
                                fclaw_bool time_interp,
                                int igrid,
                                fclaw2d_transform_data_t* transform_data);

void fclaw2d_patch_copy_face(fclaw2d_global_t* glob,
                             fclaw2d_patch_t *this_patch,
                             fclaw2d_patch_t *neighbor_patch,
                             int iface,
                             int time_interp,
                             fclaw2d_transform_data_t *transform_data);

void fclaw2d_patch_copy_corner(fclaw2d_global_t* glob,
                               fclaw2d_patch_t *this_patch,
                               fclaw2d_patch_t *corner_patch,
                               int icorner,
                               int time_interp,
                               fclaw2d_transform_data_t *transform_data);

void fclaw2d_patch_average_corner(fclaw2d_global_t* glob,
                                  fclaw2d_patch_t *coarse_patch,
                                  fclaw2d_patch_t *fine_patch,
                                  int coarse_corner,
                                  int refratio,
                                  fclaw_bool time_interp,
                                  fclaw2d_transform_data_t* transform_data);

void fclaw2d_patch_interpolate_corner(fclaw2d_global_t* glob,
                                      fclaw2d_patch_t* coarse_patch,
                                      fclaw2d_patch_t* fine_patch,
                                      int coarse_corner,
                                      int refratio,
                                      fclaw_bool time_interp,
                                      fclaw2d_transform_data_t* transform_data);

int fclaw2d_patch_tag4refinement(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *this_patch,
                                 int blockno, int patchno,
                                 int initflag);

int fclaw2d_patch_tag4coarsening(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *fine_patches,
                                 int blockno,
                                 int patchno);

void fclaw2d_patch_interpolate2fine(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t* coarse_patch,
                                    fclaw2d_patch_t* fine_patches,
                                    int this_blockno, int coarse_patchno,
                                    int fine0_patchno);

void fclaw2d_patch_average2coarse(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *fine_patches,
                                  fclaw2d_patch_t *coarse_patch,
                                  int blockno, int fine0_patchno,
                                  int coarse_patchno);

void fclaw2d_patch_setup_timeinterp(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    double alpha);

#if 0
void fclaw2d_patch_compute_diagnostics(fclaw2d_domain_t* domain, fclaw2d_patch_t* this_patch,
                                       int this_block_idx, int this_patch_idx,
                                       void* local_accumulator);
#endif 
void fclaw2d_patch_write_file(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              int iframe,int patch_num,
                              int level);

double fclaw2d_patch_single_step_update(fclaw2d_global_t *glob,
                                        fclaw2d_patch_t *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        double t,
                                        double dt);



typedef void* (*fclaw2d_patch_new_t)();

typedef void (*fclaw2d_patch_delete_t)(void *user_patch);

typedef void (*fclaw2d_patch_delete_ghost_t)(void *user_patch);

typedef void (*fclaw2d_patch_setup_t)(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx);

typedef void (*fclaw2d_patch_setup_ghost_t)(fclaw2d_global_t *glob,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx);

typedef void (*fclaw2d_patch_initialize_t)(fclaw2d_global_t *glob,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx);

typedef void (*fclaw2d_patch_physical_bc_t)(fclaw2d_global_t *glob,
                                            fclaw2d_patch_t *this_patch,
                                            int this_block_idx,
                                            int this_patch_idx,
                                            double t,
                                            double dt,
                                            fclaw_bool *intersects_bc,
                                            fclaw_bool time_interp);

typedef double (*fclaw2d_patch_single_step_update_t)(fclaw2d_global_t *glob,
                                                     fclaw2d_patch_t *this_patch,
                                                     int this_block_idx,
                                                     int this_patch_idx,
                                                     double t,
                                                     double dt);

typedef void (*fclaw2d_patch_copy_face_t)(fclaw2d_global_t* glob,
                                          fclaw2d_patch_t *this_patch,
                                          fclaw2d_patch_t *neighbor_patch,
                                          int iface,
                                          int time_interp,
                                          fclaw2d_transform_data_t *transform_data);

typedef void (*fclaw2d_patch_average_face_t)(fclaw2d_global_t* glob,
                                             fclaw2d_patch_t *coarse_patch,
                                             fclaw2d_patch_t *fine_patch,
                                             int idir,
                                             int iface_coarse,
                                             int p4est_refineFactor,
                                             int refratio,
                                             fclaw_bool time_interp,
                                             int igrid,
                                             fclaw2d_transform_data_t* transform_data);

typedef void (*fclaw2d_patch_interpolate_face_t)(fclaw2d_global_t* glob,
                                                 fclaw2d_patch_t *coarse_patch,
                                                 fclaw2d_patch_t *fine_patch,
                                                 int idir,
                                                 int iside,
                                                 int p4est_refineFactor,
                                                 int refratio,
                                                 fclaw_bool a_time_interp,
                                                 int igrid,
                                                 fclaw2d_transform_data_t* transform_data);

typedef void (*fclaw2d_patch_copy_corner_t)(fclaw2d_global_t* glob,
                                            fclaw2d_patch_t *this_patch,
                                            fclaw2d_patch_t *corner_patch,
                                            int icorner,
                                            int time_interp,
                                            fclaw2d_transform_data_t *transform_data);

typedef void (*fclaw2d_patch_average_corner_t)(fclaw2d_global_t* glob,
                                               fclaw2d_patch_t *coarse_patch,
                                               fclaw2d_patch_t *fine_patch,
                                               int coarse_corner,
                                               int refratio,
                                               fclaw_bool time_interp,
                                               fclaw2d_transform_data_t* transform_data);

typedef void (*fclaw2d_patch_interpolate_corner_t)(fclaw2d_global_t* glob,
                                                   fclaw2d_patch_t* coarse_patch,
                                                   fclaw2d_patch_t* fine_patch,
                                                   int coarse_corner,
                                                   int refratio,
                                                   fclaw_bool a_time_interp,
                                                   fclaw2d_transform_data_t* transform_data);

typedef int (*fclaw2d_patch_tag4refinement_t)(fclaw2d_global_t *glob,
                                              fclaw2d_patch_t *this_patch,
                                              int this_block_idx, int this_patch_idx,
                                              int initflag);

typedef int (*fclaw2d_patch_tag4coarsening_t)(fclaw2d_global_t *glob,
                                               fclaw2d_patch_t *this_patch,
                                               int this_blockno,
                                               int this_patchno);

typedef void (*fclaw2d_patch_interpolate2fine_t)(fclaw2d_global_t *glob,
                                                 fclaw2d_patch_t *coarse_patch,
                                                 fclaw2d_patch_t* fine_patches,
                                                 int this_blockno, int coarse_patchno,
                                                 int fine_patchno);

typedef void (*fclaw2d_patch_average2coarse_t)(fclaw2d_global_t *glob,
                                               fclaw2d_patch_t *fine_siblings,
                                               fclaw2d_patch_t *coarse_patch,
                                               int blockno, int fine_patchno,
                                               int coarse_patchno);

typedef void (*fclaw2d_patch_write_header_t)(fclaw2d_global_t* glob,
                                             int iframe);

typedef void (*fclaw2d_patch_write_file_t)(fclaw2d_global_t *glob,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx,
                                           int iframe,int patch_num,
                                           int level);

typedef void (*fclaw2d_patch_ghost_pack_t)(fclaw2d_global_t *glob,
                                           fclaw2d_patch_t *this_patch,
                                           double *patch_data,
                                           int time_interp);

typedef void (*fclaw2d_patch_ghost_unpack_t)(fclaw2d_global_t *glob,
                                             fclaw2d_patch_t* this_patch,
                                             int this_block_idx, int this_patch_idx,
                                             double *qdata, fclaw_bool time_interp);

typedef void (*fclaw2d_patch_build_ghost_t)(fclaw2d_global_t *glob,
                                            fclaw2d_patch_t *this_patch,
                                            int blockno,
                                            int patchno,
                                            void *user);

typedef size_t (*fclaw2d_patch_ghost_packsize_t)(fclaw2d_global_t* glob);

typedef void (*fclaw2d_patch_local_ghost_alloc_t)(fclaw2d_global_t* glob,
                                                 fclaw2d_patch_t* this_patch,
                                                 void** q);

typedef void (*fclaw2d_patch_local_ghost_free_t)(fclaw2d_global_t* glob,
                                                 void **q);

typedef void (*fclaw2d_patch_partition_pack_t)(fclaw2d_global_t *glob,
                                               fclaw2d_patch_t *this_patch,
                                               int this_block_idx,
                                               int this_patch_idx,
                                               void *user);

typedef void (*fclaw2d_patch_partition_unpack_t)(fclaw2d_global_t *glob,
                                                 fclaw2d_domain_t *new_domain,
                                                 fclaw2d_patch_t *this_patch,
                                                 int this_block_idx,
                                                 int this_patch_idx,
                                                 void *user);

typedef size_t (*fclaw2d_patch_partition_packsize_t)(fclaw2d_global_t* glob);

typedef void (*fclaw2d_patch_build_t)(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno,
                                      void *user);

typedef void (*fclaw2d_patch_build_from_fine_t)(fclaw2d_global_t *glob,
                                                fclaw2d_patch_t *fine_patches,
                                                fclaw2d_patch_t *coarse_patch,
                                                int blockno,
                                                int coarse_patchno,
                                                int fine0_patchno,
                                                fclaw2d_build_mode_t build_mode);

typedef void (*fclaw2d_patch_setup_timeinterp_t)(fclaw2d_global_t *glob,
                                                 fclaw2d_patch_t *this_patch,
                                                 double alpha);


typedef void (*fclaw2d_patch_restore_step_t)(fclaw2d_global_t *glob,
                                             fclaw2d_patch_t* this_patch);

typedef void (*fclaw2d_patch_save_step_t)(fclaw2d_global_t *glob,
                                          fclaw2d_patch_t* this_patch);


typedef struct fclaw2d_patch_vtable
{
    fclaw2d_patch_new_t                patch_new;
    fclaw2d_patch_delete_t             patch_delete;
    fclaw2d_patch_delete_ghost_t       delete_ghost;
    fclaw2d_patch_setup_t              setup;
    fclaw2d_patch_setup_ghost_t        setup_ghost;
    fclaw2d_patch_initialize_t         initialize;
    fclaw2d_patch_physical_bc_t        physical_bc;
    fclaw2d_patch_single_step_update_t single_step_update;
    fclaw2d_patch_build_t              build;
    fclaw2d_patch_build_from_fine_t    build_from_fine;
    fclaw2d_patch_restore_step_t       restore_step;
    fclaw2d_patch_save_step_t          save_step;

    /* regridding functions */
    fclaw2d_patch_tag4refinement_t     tag4refinement;
    fclaw2d_patch_tag4coarsening_t     tag4coarsening;
    fclaw2d_patch_average2coarse_t     average2coarse;
    fclaw2d_patch_interpolate2fine_t   interpolate2fine;

    /* ghost filling functions */
    fclaw2d_patch_copy_face_t           copy_face;
    fclaw2d_patch_average_face_t        average_face;
    fclaw2d_patch_interpolate_face_t    interpolate_face;

    fclaw2d_patch_copy_corner_t         copy_corner;
    fclaw2d_patch_average_corner_t      average_corner;
    fclaw2d_patch_interpolate_corner_t  interpolate_corner;

    /* output functions */
    fclaw2d_patch_write_header_t       write_header;
    fclaw2d_patch_write_file_t         write_file;

    /* Time interpolation */
    fclaw2d_patch_setup_timeinterp_t   setup_timeinterp;

    /* Ghost packing functions (for parallel use) */
    fclaw2d_patch_ghost_pack_t         ghost_pack;
    fclaw2d_patch_ghost_unpack_t       ghost_unpack;
    fclaw2d_patch_build_ghost_t        build_ghost;
    fclaw2d_patch_ghost_packsize_t     ghost_packsize;
    fclaw2d_patch_local_ghost_alloc_t  local_ghost_alloc;
    fclaw2d_patch_local_ghost_free_t   local_ghost_free;

    /* partitioning */
    fclaw2d_patch_partition_pack_t       partition_pack;
    fclaw2d_patch_partition_unpack_t     partition_unpack;
    fclaw2d_patch_partition_packsize_t   partition_packsize;

    int defaults_set;

} fclaw2d_patch_vtable_t;

void fclaw2d_set_patch_vtable(fclaw2d_patch_vtable_t user_vt);
fclaw2d_patch_vtable_t* fclaw2d_patch_vt();


/* -----------------------------------------------------
   Ghost exchange
   ---------------------------------------------------- */
#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
