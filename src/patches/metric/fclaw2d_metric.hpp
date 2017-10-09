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

#ifndef FCLAW2D_METRIC_PATCH_HPP
#define FCLAW2D_METRIC_PATCH_HPP

#include <fclaw2d_farraybox.hpp>

#include <fclaw2d_patch.h>                /* Needed to get enum for build modes */

#include "fclaw2d_metric_default_fort.h"  /* Needed for fort typdefs in vtable */

typedef struct fclaw2d_metric_vtable fclaw2d_metric_vtable_t;

class fclaw2d_clawpatch_t;
struct fclaw2d_global;
struct fclaw2d_patch;

class fclaw2d_metric_patch_t
{
public :

    int mx;           
    int my;           
    int mbc;
    int blockno;

    double dx;
    double dy;
    double xlower;
    double ylower;
    double xupper;
    double yupper;
    
    FArrayBox xp;
    FArrayBox yp;
    FArrayBox zp;

    FArrayBox xd;
    FArrayBox yd;
    FArrayBox zd;

    FArrayBox xface_normals;
    FArrayBox yface_normals;
    FArrayBox xface_tangents;
    FArrayBox yface_tangents;
    FArrayBox surf_normals;
    FArrayBox edge_lengths;

    FArrayBox area;
    FArrayBox curvature;  // ???
};

/* --------------------------- Metric routines (typedefs) ----------------------------- */

typedef void (*fclaw2d_metric_compute_mesh_t)(struct fclaw2d_global *glob,
                                              struct fclaw2d_patch *this_patch,
                                              int blockno,
                                              int patchno);

typedef void (*fclaw2d_metric_compute_area_t)(struct fclaw2d_global *glob,
                                              struct fclaw2d_patch *patch,
                                              int blockno,
                                              int patchno);

typedef void (*fclaw2d_metric_compute_area_ghost_t)(struct fclaw2d_global *glob,
                                                    struct fclaw2d_patch *patch,
                                                    int blockno,
                                                    int patchno);

typedef void (*fclaw2d_metric_compute_normals_t)(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *this_patch,
                                                 int blockno,
                                                 int patchno);


/* ------------------------------ Metric routines ------------------------------------- */

fclaw2d_metric_patch_t* fclaw2d_metric_patch_new();

void fclaw2d_metric_patch_delete(fclaw2d_metric_patch_t *patchmp);


void fclaw2d_metric_patch_define(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch *this_patch,
                                 int blockno, int patchno,
                                 fclaw2d_build_mode_t build_mode);


void fclaw2d_metric_patch_setup(struct fclaw2d_global* glob,
                                struct fclaw2d_patch* this_patch,
                                int blockno,
                                int patchno);

void fclaw2d_metric_patch_setup_from_fine(struct fclaw2d_global *glob,
                                          struct fclaw2d_patch *fine_patches,
                                          struct fclaw2d_patch *coarse_coarse,
                                          int blockno,
                                          int coarse_patchno,
                                          int fine0_patchno);

void fclaw2d_metric_patch_compute_area(struct fclaw2d_global *glob,
                                       struct fclaw2d_patch* this_patch,
                                       int blockno, int patchno);

/* --------------------------------- Access functions --------------------------------- */

void fclaw2d_metric_patch_grid_data(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* this_patch,
                                    int* mx, int* my, int* mbc,
                                    double* xlower, double* ylower,
                                    double* dx, double* dy);


void fclaw2d_metric_patch_mesh_data(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* this_patch,
                                    double **xp, double **yp, double **zp,
                                    double **xd, double **yd, double **zd,
                                    double **area);


void fclaw2d_metric_patch_mesh_data2(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* this_patch,
                                     double **xnormals, double **ynormals,
                                     double **xtangents, double **ytangents,
                                     double **surfnormals,
                                     double **edgelengths, double **curvature);


double* fclaw2d_metric_patch_get_area(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch* this_patch);


/* ---------------------------- Metric default (virtualized) -------------------------- */

void fclaw2d_metric_compute_area_default(struct fclaw2d_global *glob,
                                         struct fclaw2d_patch* this_patch,
                                         int blockno, int patchno);


void fclaw2d_metric_compute_area_ghost_default(struct fclaw2d_global* glob,
                                               struct fclaw2d_patch* this_patch,
                                               int blockno,
                                               int patchno);

void fclaw2d_metric_compute_mesh_default(struct fclaw2d_global *glob,
                                         struct fclaw2d_patch *this_patch,
                                         int blockno,
                                         int patchno);


void fclaw2d_metric_compute_normals_default(struct fclaw2d_global *glob,
                                            struct fclaw2d_patch *this_patch,
                                            int blockno,
                                            int patchno);


/* ------------------------------------ Virtual table --------------------------------- */

fclaw2d_metric_vtable_t* fclaw2d_metric_vt();


struct fclaw2d_metric_vtable
{
    /* Building patches, including functions to create metric terms */
    fclaw2d_metric_compute_mesh_t        compute_mesh;    /* wrapper */
    fclaw2d_metric_compute_area_t        compute_area;  /* wrapper */
    fclaw2d_metric_compute_area_ghost_t  compute_area_ghost;
    fclaw2d_metric_compute_normals_t     compute_normals;  /* wrapper */

    /* Fortran files */
    fclaw2d_fort_compute_mesh_t          fort_compute_mesh;
    fclaw2d_fort_compute_normals_t       fort_compute_normals;
    fclaw2d_fort_compute_tangents_t      fort_compute_tangents;
    fclaw2d_fort_compute_surf_normals_t  fort_compute_surf_normals;

    int is_set;
};

void fclaw2d_metric_vtable_initialize();




#endif /* !FCLAW2D_METRIC_PATCH_HPP */
