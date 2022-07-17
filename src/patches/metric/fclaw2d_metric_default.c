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

#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <fclaw_math.h>

#include "fclaw2d_metric.h"
#include "fclaw2d_metric_default_fort.h"


void fclaw2d_metric_compute_mesh_default(fclaw2d_global_t *glob,
                                         fclaw2d_patch_t* patch,
                                         int blockno,
                                         int patchno)
{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt(glob);
    FCLAW_ASSERT(metric_vt->fort_compute_mesh);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_metric_patch_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *area;
    fclaw2d_metric_patch_mesh_data(glob,patch,
                                         &xp,&yp,&zp,&xd,&yd,&zd,&area);

    /* Compute centers and corners of mesh cell */
    metric_vt->fort_compute_mesh(&mx,&my,&mbc,&xlower,&ylower,&dx,&dy,
                                 &blockno,xp,yp,zp,xd,yd,zd);
}


void fclaw2d_metric_compute_tensors_default(fclaw2d_global_t *glob,
                                            fclaw2d_patch_t *patch,
                                            int blockno,
                                            int patchno)
{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt(glob);


    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_metric_patch_grid_data(glob,patch,&mx,&my,&mbc,
                                   &xlower,&ylower,&dx,&dy);

    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *area;
    fclaw2d_metric_patch_mesh_data(glob,patch,
                                   &xp,&yp,&zp,&xd,&yd,&zd,&area);

    double *xnormals, *ynormals;
    double *xtangents, *ytangents;
    double *edgelengths;
    double *surfnormals, *curvature;
    fclaw2d_metric_patch_mesh_data2(glob,patch,
                                    &xnormals,&ynormals,
                                    &xtangents,&ytangents,
                                    &surfnormals,&edgelengths,
                                    &curvature);

    /* The user could set these to NULL to avoid doing these computations ... */

    if (metric_vt->fort_compute_normals != NULL)
    {
        metric_vt->fort_compute_normals(&mx,&my,&mbc,xp,yp,zp,xd,yd,zd,
                                        xnormals,ynormals);
    }

    if (metric_vt->fort_compute_tangents != NULL)
    {
        metric_vt->fort_compute_tangents(&mx,&my,&mbc,xd,yd,zd,xtangents,ytangents,
                                         edgelengths);
    }

    if (metric_vt->fort_compute_surf_normals != NULL)
    {
        metric_vt->fort_compute_surf_normals(&mx,&my,&mbc,xnormals,ynormals,edgelengths,
                                             curvature, surfnormals, area);
    }
}


void fclaw2d_metric_compute_area_default(fclaw2d_global_t *glob,
                                         fclaw2d_patch_t *patch,
                                         int blockno,
                                         int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_metric_patch_grid_data(glob,patch,&mx,&my,&mbc,
                                   &xlower,&ylower,&dx,&dy);

    double *area = fclaw2d_metric_patch_get_area(glob, patch);

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    int level = patch->level;
    int maxlevel = fclaw_opt->maxlevel;
    int refratio = fclaw_opt->refratio;
    int m = pow_int(refratio,maxlevel-level);
    double *quadstore = FCLAW_ALLOC(double,3*(m+1)*(m+1));

    /* Fix this so we allocate quadstore in the Fortran routine */
    int ghost_only = 0;
    FCLAW2D_FORT_COMPUTE_AREA(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                              &blockno, area, &m, quadstore, &ghost_only);

    FCLAW_FREE(quadstore);
}

void fclaw2d_metric_compute_area_ghost_default(fclaw2d_global_t* glob,
                                               fclaw2d_patch_t* patch,
                                               int blockno,
                                               int patchno)
{

    int mx,my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_metric_patch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_metric_patch_get_area(glob, patch);

    /* Set area in ghost cells not set above */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    int level = patch->level;
    int maxlevel = fclaw_opt->maxlevel;
    int refratio = fclaw_opt->refratio;

    int m = pow_int(refratio,maxlevel-level);
    double *quadstore = FCLAW_ALLOC(double,3*(m+1)*(m+1));

    double *area;
    int ghost_only = 1;
    FCLAW2D_FORT_COMPUTE_AREA(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                              &blockno, area, &m, quadstore,
                              &ghost_only);
    FCLAW_FREE(quadstore);
}