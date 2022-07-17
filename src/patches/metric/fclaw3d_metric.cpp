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


#define REFINE_DIM 2
#define PATCH_DIM 3

#include <fclaw2d_metric.cpp>


void fclaw3d_metric_patch_scalar(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* patch,
                                 double **volume, double** facearea)
{
    fclaw3d_metric_patch_t* mp = fclaw3d_metric_get_metric_patch(glob, patch);
    *vol = mp->volume.dataPtr();
    *faceareas =  mp->face_area.dataPtr();
}


void fclaw3d_metric_patch_basis(fclaw2d_global_t* glob,
                                fclaw2d_patch_t* patch,
                                double **xrot, double **yrot, double **zrot)
{
    fclaw3d_metric_patch_t* mp = fclaw3d_metric_get_metric_patch(glob, patch);
    *xrot = mp->xrot.dataPtr();
    *yrot = mp->yrot.dataPtr();
    *zrot = mp->zrot.dataPtr();
}


void fclaw3d_metric_patch_grid_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* patch,
                                    int* mx, int* my, int* mz, 
                                    int* mbc,
                                    double* xlower, double* ylower, double zlower,
                                    double* dx, double* dy, double* dz)
{
    fclaw3d_metric_patch_t* mp = fclaw3d_metric_get_metric_patch(glob, patch);
    *mx     = mp->mx;
    *my     = mp->my;
    *mz     = mp->mz;
    *mbc    = mp->mbc;
    *xlower = mp->xlower;
    *ylower = mp->ylower;
    *zlower = mp->zlower;
    *dx     = mp->dx;
    *dy     = mp->dy;
    *dz     = mp->dz;
}


void fclaw3d_metric_patch_mesh_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* patch,
                                    double **xp, double **yp, double **zp,
                                    double **xd, double **yd, double **zd,
                                    double **volume)
{
    fclaw3d_metric_patch_t* mp = fclaw3d_metric_get_metric_patch(glob, patch);
    *xp = mp->xp.dataPtr();
    *yp = mp->yp.dataPtr();
    *zp = mp->zp.dataPtr();
    *xd = mp->xd.dataPtr();
    *yd = mp->yd.dataPtr();
    *zd = mp->zd.dataPtr();
    *volume = mp->volume.dataPtr();
}

void fclaw3d_metric_patch_mesh_data2(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* patch,
                                     double **xrot, double **yrot, 
                                     double **zrot,
                                     double **faceareas)
{
    fclaw3d_metric_patch_t* mp = fclaw3d_metric_get_metric_patch(glob, patch);
    *xrot = mp->xrot;
    *yrot = mp->yrot;
    *zrot = mp->zrot;
    *faceareas = mp->faceareas;
}


