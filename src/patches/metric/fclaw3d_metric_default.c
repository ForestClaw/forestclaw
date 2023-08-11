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

#include <fclaw_global.h>
#include <fclaw2d_options.h>
#include <fclaw_math.h>

#include "fclaw3d_metric.h"
#include "fclaw3d_metric_default_fort.h"


void fclaw3d_metric_compute_mesh_default(fclaw_global_t *glob,
                                         fclaw_patch_t* patch,
                                         int blockno,
                                         int patchno)
{
    fclaw3d_metric_vtable_t *metric_vt = fclaw3d_metric_vt(glob);

    int mx,my,mz,mbc;
    double xlower,ylower,zlower,dx,dy,dz;
    fclaw3d_metric_patch_grid_data(glob,patch,
                                   &mx,&my,&mz,&mbc,
                                   &xlower,&ylower,&zlower,
                                   &dx,&dy,&dz);
 
    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *volume, *facearea; /* Not used here */
    fclaw3d_metric_patch_mesh_data(glob,patch,
                                   &xp,&yp,&zp,&xd,&yd,&zd,&volume,&facearea);

    /* Compute centers and corners of mesh cell */
    FCLAW_ASSERT(metric_vt->fort_compute_mesh != NULL);
    metric_vt->fort_compute_mesh(&mx,&my,&mz,&mbc,&xlower,&ylower,&zlower,
                                 &dx,&dy,&dz,
                                 &blockno,xp,yp,zp,xd,yd,zd);
}


void fclaw3d_metric_compute_basis_default(fclaw_global_t *glob,
                                          fclaw_patch_t *patch,
                                          int blockno,
                                          int patchno)
{
    fclaw3d_metric_vtable_t *metric_vt = fclaw3d_metric_vt(glob);


    int mx,my,mz,mbc;
    double xlower,ylower,zlower,dx,dy,dz;
    fclaw3d_metric_patch_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                   &xlower,&ylower,&zlower, &dx,&dy, &dz);

    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *volume, *facearea;  
    fclaw3d_metric_patch_mesh_data(glob,patch,
                                   &xp,&yp,&zp,&xd,&yd,&zd,
                                   &volume,&facearea);

    /* The user could set these to NULL to avoid doing these computations ... */

    if (metric_vt->fort_compute_basis != NULL)
    {
        double *xrot, *yrot, *zrot;
        fclaw3d_metric_patch_basis(glob,patch,&xrot,&yrot,&zrot);

        int ghost_only = 0;
        metric_vt->fort_compute_basis(&mx,&my,&mz,&mbc,xd,yd,zd,
                                       xrot,yrot,zrot,&ghost_only);
    }
}


void fclaw3d_metric_compute_volume_default(fclaw_global_t *glob,
                                         fclaw_patch_t *patch,
                                         int blockno,
                                         int patchno)
{
    int mx,my,mz,mbc;
    double xlower,ylower,zlower,dx,dy,dz;
    fclaw3d_metric_patch_grid_data(glob,patch,&mx,&my,&mz,&mbc,
                                   &xlower,&ylower,&zlower,&dx,&dy,&dz);

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    int level = patch->level;
    int maxlevel = fclaw_opt->maxlevel;
    int refratio = 2;

    /* Set up a local fine grid for each cell;   compute volumes on fine
       grid sub-cells;  coarse grid volume is then sum of fine grid volumes.
       This is needed for geometric consistency 
    */
    int m = pow_int(refratio,maxlevel-level);
    double *hexstore = FCLAW_ALLOC(double,3*(m+1)*(m+1)*(m+1));

    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *volume, *facearea;
    fclaw3d_metric_patch_mesh_data(glob,patch,
                                   &xp,&yp,&zp,&xd,&yd,&zd,
                                   &volume,&facearea);

    int xd_size = fclaw3d_metric_patch_nodes_size(glob,patch);
    if (xd_size == 0)
    {
        printf("xd size is 0\n");
        exit(0);
    }
    /* Computes volume and face areas */
    int ghost_only = 0;

    FCLAW_ASSERT(xd != NULL);  
    FCLAW_ASSERT(facearea != NULL);
    FCLAW3D_METRIC_FORT_COMPUTE_VOLUME(&mx, &my, &mz, &mbc, &blockno,
                                       &dx, &dy, &dz,
                                       &xlower, &ylower,&zlower,
                                       xd,yd,zd,
                                       volume, facearea,
                                       &m, hexstore,
                                       &ghost_only);

    FCLAW_FREE(hexstore);
}

void fclaw3d_metric_compute_volume_ghost_default(fclaw_global_t* glob,
                                                 fclaw_patch_t* patch,
                                                 int blockno,
                                                 int patchno)
{
    /* This is used for filling coarse grid ghost cells when the coarse grid
       has been constructed by averaging solution and metric terms from a fine 
       grid.  In this case, we only average internal cell values;  ghost values
       do not get averaged.
    */

    int mx,my, mz, mbc;
    double xlower,ylower,zlower, dx,dy, dz;
    fclaw3d_metric_patch_grid_data(glob,patch,&mx,&my,&mz,&mbc,
                                   &xlower,&ylower,&zlower,
                                   &dx, &dy, &dz);


    /* Set area in ghost cells not set above */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    int level = patch->level;
    int maxlevel = fclaw_opt->maxlevel;
    int refratio = fclaw_opt->refratio;

    int m = pow_int(refratio,maxlevel-level);
    double *hexstore = FCLAW_ALLOC(double,3*(m+1)*(m+1)*(m+1));

    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *volume, *faceareas;
    fclaw3d_metric_patch_mesh_data(glob,patch,
                                   &xp,&yp,&zp,&xd,&yd,&zd,
                                   &volume,&faceareas);

    int ghost_only = 1;
    FCLAW3D_METRIC_FORT_COMPUTE_VOLUME(&mx, &my, &mz, &mbc, &blockno,
                                &dx, &dy, &dz, 
                                &xlower, &ylower,&zlower,
                                xd,yd,zd,
                                volume, faceareas, &m, hexstore,
                                &ghost_only);
    FCLAW_FREE(hexstore);
}