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

/* Needed as long as metric type defs are defined here */
#include <fclaw2d_metric_default.h>  
#include <fclaw2d_metric_default_fort.h>


#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

typedef void (*fclaw2d_vtable_initialize_t)();

typedef void (*fclaw2d_problem_setup_t)(struct fclaw2d_global *glob);

typedef void (*fclaw2d_output_frame_t)(struct fclaw2d_global * glob, int iframe);


typedef void (*fclaw2d_metric_setup_mesh_t)(struct fclaw2d_global *glob,
                                            struct fclaw2d_patch *this_patch,
                                            int blockno,
                                            int patchno);

typedef void (*fclaw2d_metric_compute_area_t)(struct fclaw2d_global *glob,
                                              struct fclaw2d_patch* this_patch,
                                              int blockno,
                                              int patchno);

typedef void (*fclaw2d_metric_area_set_ghost_t)(struct fclaw2d_global *glob,
                                                struct fclaw2d_patch* this_patch,
                                                int blockno,
                                                int patchno);

typedef void (*fclaw2d_metric_compute_normals_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *this_patch,
                                                 int blockno,
                                                 int patchno);

typedef void (*fclaw2d_after_regrid_t)(struct fclaw2d_global *glob);

typedef struct fclaw2d_vtable
{

    fclaw2d_vtable_initialize_t          vtable_init;

    fclaw2d_problem_setup_t              problem_setup;

    /* regridding functions */
    fclaw2d_after_regrid_t               after_regrid;

    /* Output functions */
    fclaw2d_output_frame_t               output_frame;

    /* Building patches, including functions to create metric terms */
    fclaw2d_metric_setup_mesh_t          metric_setup_mesh;    /* wrapper */
    fclaw2d_fort_setup_mesh_t            fort_setup_mesh;

    fclaw2d_metric_compute_area_t        metric_compute_area;  /* wrapper */
    fclaw2d_metric_area_set_ghost_t      metric_area_set_ghost;

    fclaw2d_metric_compute_normals_t     metric_compute_normals;  /* wrapper */
    fclaw2d_fort_compute_normals_t       fort_compute_normals;
    fclaw2d_fort_compute_tangents_t      fort_compute_tangents;
    fclaw2d_fort_compute_surf_normals_t  fort_compute_surf_normals;

    int is_set;


} fclaw2d_vtable_t;


fclaw2d_vtable_t* fclaw2d_vt();
void fclaw2d_vtable_initialize();


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
