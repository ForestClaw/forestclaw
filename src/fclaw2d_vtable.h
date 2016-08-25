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

#ifndef FCLAW2D_VTABLE_H
#define FCLAW2D_VTABLE_H

#include <forestclaw2d.h>
#include <fclaw_base.h>
#include <fclaw2d_defs.h>

#include <fclaw2d_output.h>
#include <fclaw2d_output_ascii.h>
#include <fclaw2d_regrid_default.h>
#include <fclaw2d_metric_default.h>
#include <fclaw2d_diagnostics_default.h>
#include <fclaw2d_transform.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


typedef void (*fclaw2d_problem_setup_t)(fclaw2d_domain_t *domain);

typedef void (*fclaw2d_patch_setup_t)(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx);

typedef void (*fclaw2d_patch_initialize_t)(fclaw2d_domain_t *domain,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx);

typedef void (*fclaw2d_patch_physical_bc_t)(fclaw2d_domain_t *domain,
                                            fclaw2d_patch_t *this_patch,
                                            int this_block_idx,
                                            int this_patch_idx,
                                            double t,
                                            double dt,
                                            fclaw_bool *intersects_bc,
                                            fclaw_bool time_interp);

typedef double (*fclaw2d_patch_single_step_update_t)(fclaw2d_domain_t *domain,
                                                     fclaw2d_patch_t *this_patch,
                                                     int this_block_idx,
                                                     int this_patch_idx,
                                                     double t,
                                                     double dt);

typedef void (*fclaw2d_metric_setup_mesh_t)(fclaw2d_domain_t *domain,
                                            fclaw2d_patch_t *this_patch,
                                            int blockno,
                                            int patchno);

typedef void (*fclaw2d_metric_compute_area_t)(fclaw2d_domain_t *domain,
                                              fclaw2d_patch_t* this_patch,
                                              int blockno,
                                              int patchno);

typedef void (*fclaw2d_metric_area_set_ghost_t)(fclaw2d_domain_t *domain,
                                                fclaw2d_patch_t* this_patch,
                                                int blockno,
                                                int patchno);

typedef void (*fclaw2d_metric_compute_normals_t)(fclaw2d_domain_t *domain,
                                                 fclaw2d_patch_t *this_patch,
                                                 int blockno,
                                                 int patchno);

typedef int (*fclaw2d_regrid_tag4refinement_t)(fclaw2d_domain_t *domain,
                                              fclaw2d_patch_t *this_patch,
                                              int this_block_idx, int this_patch_idx,
                                              int initflag);

typedef int (*fclaw2d_regrid_tag4coarsening_t)(fclaw2d_domain_t *domain,
                                               fclaw2d_patch_t *this_patch,
                                               int this_blockno,
                                               int this_patchno);

typedef void (*fclaw2d_regrid_interpolate2fine_t)(fclaw2d_domain_t* domain,
                                                 fclaw2d_patch_t *coarse_patch,
                                                 fclaw2d_patch_t* fine_patches,
                                                 int this_blockno, int coarse_patchno,
                                                 int fine_patchno);

typedef void (*fclaw2d_regrid_average2coarse_t)(fclaw2d_domain_t *domain,
                                               fclaw2d_patch_t *fine_siblings,
                                               fclaw2d_patch_t *coarse_patch,
                                               int blockno, int fine_patchno,
                                               int coarse_patchno);

typedef void (*fclaw2d_patch_write_header_t)(fclaw2d_domain_t* domain,
                                             int iframe);

typedef void (*fclaw2d_patch_write_file_t)(fclaw2d_domain_t *domain,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx,
                                           int iframe,int patch_num,
                                           int level);



typedef void (*fclaw2d_run_user_diagnostics_t)(fclaw2d_domain_t *domain, const double t);
typedef void (*fclaw2d_diagnostics_compute_error_t)(fclaw2d_domain_t *domain,
                                                    fclaw2d_patch_t *this_patch,
                                                    int this_block_idx,
                                                    int this_patch_idx,
                                                    double *error);

typedef void (*fclaw2d_fort_compute_error_t)(int* blockno, int *mx, int *my, int *mbc,
                                             int *meqn,
                                             double *dx, double *dy, double *xlower,
                                             double *ylower, double *t, double q[],
                                             double error[]);

typedef void (*fclaw2d_patch_copy_face_t)(fclaw2d_domain_t *domain,
                                                fclaw2d_patch_t *this_patch,
                                                fclaw2d_patch_t *neighbor_patch,
                                                int iface,
                                                int time_interp,
                                                fclaw2d_transform_data_t *transform_data);

typedef void (*fclaw2d_patch_copy_corner_t)(fclaw2d_domain_t *domain,
                                                  fclaw2d_patch_t *this_patch,
                                                  fclaw2d_patch_t *corner_patch,
                                                  int icorner,
                                                  int time_interp,
                                                  fclaw2d_transform_data_t *transform_data);

typedef void (*fclaw2d_fort_copy_face_t)(int* mx, int* my,int* mbc, int* meqn,
                                               double qthis[],double qneighbor[],
                                               int* iface,
                                               fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw2d_fort_copy_corner_t)(int* mx, int* my, int* mbc, int* meqn,
                                                 double this_q[],double neighbor_q[],
                                                 int* icorner,
                                                 fclaw2d_transform_data_t** transform_cptr);


typedef struct fclaw2d_vtable
{
    fclaw2d_problem_setup_t            problem_setup;

    fclaw2d_patch_setup_t              patch_setup;
    fclaw2d_patch_initialize_t         patch_initialize;
    fclaw2d_patch_physical_bc_t        patch_physical_bc;
    fclaw2d_patch_single_step_update_t patch_single_step_update;

    /* Building patches, including functions to create metric terms */
    fclaw2d_metric_setup_mesh_t        metric_setup_mesh;    /* wrapper */
    fclaw2d_fort_setup_mesh_t          fort_setup_mesh;

    fclaw2d_metric_compute_area_t      metric_compute_area;  /* wrapper */
    fclaw2d_metric_area_set_ghost_t    metric_area_set_ghost;

    fclaw2d_metric_compute_normals_t    metric_compute_normals;  /* wrapper */
    fclaw2d_fort_compute_normals_t      fort_compute_normals;
    fclaw2d_fort_compute_tangents_t     fort_compute_tangents;
    fclaw2d_fort_compute_surf_normals_t fort_compute_surf_normals;


    /* regridding functions */
    fclaw2d_regrid_average2coarse_t    regrid_average2coarse;
    fclaw2d_fort_average2coarse_t      fort_average2coarse;

    fclaw2d_regrid_interpolate2fine_t  regrid_interpolate2fine;
    fclaw2d_fort_interpolate2fine_t    fort_interpolate2fine;

    fclaw2d_regrid_tag4refinement_t    regrid_tag4refinement;
    fclaw2d_fort_tag4refinement_t      fort_tag4refinement;

    fclaw2d_regrid_tag4coarsening_t    regrid_tag4coarsening;
    fclaw2d_fort_tag4coarsening_t      fort_tag4coarsening;

    /* output functions */
    fclaw2d_patch_write_header_t       write_header;
    fclaw2d_fort_write_header_t        fort_write_header;

    fclaw2d_patch_write_file_t         patch_write_file;
    fclaw2d_fort_write_file_t          fort_write_file;

    /* diagnostic functions */
    fclaw2d_run_user_diagnostics_t       run_user_diagnostics;
    fclaw2d_diagnostics_compute_error_t  compute_patch_error;
    fclaw2d_fort_compute_error_t         fort_compute_patch_error;

    /* ghost filling functions */
    fclaw2d_patch_copy_face_t    copy_face;
    fclaw2d_fort_copy_face_t     fort_copy_face;

    fclaw2d_patch_copy_corner_t  copy_corner;
    fclaw2d_fort_copy_corner_t   fort_copy_corner;


} fclaw2d_vtable_t;


void
    fclaw2d_init_vtable(fclaw2d_vtable_t *vt);

void
fclaw2d_set_vtable(fclaw2d_domain_t* domain, fclaw2d_vtable_t *vt);

fclaw2d_vtable_t
fclaw2d_get_vtable(fclaw2d_domain_t *domain);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
