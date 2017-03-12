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

#include <fclaw2d_metric_default.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_transform.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


typedef void (*fclaw2d_problem_setup_t)(fclaw2d_domain_t *domain);

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

typedef void (*fclaw2d_after_regrid_t)(fclaw2d_domain_t *domain);


/* Diagnostic information */
typedef void (*fclaw2d_diagnostics_initialize_t)(fclaw2d_domain_t *domain,
                                                 void** acc);

typedef void (*fclaw2d_diagnostics_compute_t)(fclaw2d_domain_t *domain,
                                              void* acc);

typedef void (*fclaw2d_diagnostics_gather_t)(fclaw2d_domain_t *domain,
                                             void* acc,
                                             int init_flag);

typedef void (*fclaw2d_diagnostics_reset_t)(fclaw2d_domain_t *domain,
                                            void* acc);

typedef void (*fclaw2d_diagnostics_finalize_t)(fclaw2d_domain_t *domain,
                                               void** acc);

typedef struct fclaw2d_vtable
{
    fclaw2d_problem_setup_t              problem_setup;

    /* regridding functions */
    fclaw2d_after_regrid_t               after_regrid;

    /* patch diagnostic functions (error, conservation, area, etc) */
    fclaw2d_diagnostics_initialize_t     patch_init_diagnostics;
    fclaw2d_diagnostics_compute_t        patch_compute_diagnostics;
    fclaw2d_diagnostics_gather_t         patch_gather_diagnostics;
    fclaw2d_diagnostics_reset_t          patch_reset_diagnostics;
    fclaw2d_diagnostics_finalize_t       patch_finalize_diagnostics;

    /* solver diagnostic functions (gauges, fgmax, and so on) */
    fclaw2d_diagnostics_initialize_t     solver_init_diagnostics;
    fclaw2d_diagnostics_compute_t        solver_compute_diagnostics;
    fclaw2d_diagnostics_gather_t         solver_gather_diagnostics;
    fclaw2d_diagnostics_reset_t          solver_reset_diagnostics;
    fclaw2d_diagnostics_finalize_t       solver_finalize_diagnostics;

    /* user defined diagnostics */
    fclaw2d_diagnostics_initialize_t     user_init_diagnostics;
    fclaw2d_diagnostics_compute_t        user_compute_diagnostics;
    fclaw2d_diagnostics_gather_t         user_gather_diagnostics;
    fclaw2d_diagnostics_reset_t          user_reset_diagnostics;
    fclaw2d_diagnostics_finalize_t       user_finalize_diagnostics;

    /* Building patches, including functions to create metric terms */
    fclaw2d_metric_setup_mesh_t          metric_setup_mesh;    /* wrapper */
    fclaw2d_fort_setup_mesh_t            fort_setup_mesh;

    fclaw2d_metric_compute_area_t        metric_compute_area;  /* wrapper */
    fclaw2d_metric_area_set_ghost_t      metric_area_set_ghost;

    fclaw2d_metric_compute_normals_t     metric_compute_normals;  /* wrapper */
    fclaw2d_fort_compute_normals_t       fort_compute_normals;
    fclaw2d_fort_compute_tangents_t      fort_compute_tangents;
    fclaw2d_fort_compute_surf_normals_t  fort_compute_surf_normals;


} fclaw2d_vtable_t;


fclaw2d_vtable_t* fclaw2d_vt();

void fclaw2d_init_vtable();

void fclaw2d_set_vtable();

// fclaw2d_vtable_t fclaw2d_get_vtable(fclaw2d_domain_t *domain);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
