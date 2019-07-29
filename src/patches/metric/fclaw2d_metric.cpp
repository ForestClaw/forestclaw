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

#include "fclaw2d_metric.h"
#include "fclaw2d_metric.hpp"

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>  

static fclaw2d_metric_vtable_t s_metric_vt;


static
fclaw2d_metric_patch_t* get_metric_patch(fclaw2d_patch_t *this_patch)
{
    return (fclaw2d_metric_patch_t*) fclaw2d_patch_metric_patch(this_patch);
}


/* ----------------------------- Creating/deleting patches ---------------------------- */

fclaw2d_metric_patch_t* fclaw2d_metric_patch_new()
{
    fclaw2d_metric_patch_t *mp = new fclaw2d_metric_patch_t;
    return mp;
}

void fclaw2d_metric_patch_delete(fclaw2d_metric_patch_t **mp)
{
    FCLAW_ASSERT(mp != NULL);
    delete *mp;
    *mp = NULL;
}


void fclaw2d_metric_patch_define(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t *this_patch,
                                 int mx, int my, int mbc, 
                                 double dx, double dy, 
                                 double xlower, double ylower,
                                 double xupper, double yupper,
                                 int blockno, int patchno,
                                 fclaw2d_build_mode_t build_mode)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);

    mp->mx   = mx;
    mp->my   = my;
    mp->mbc  = mbc;
    mp->blockno = blockno;

    mp->dx = dx;
    mp->dy = dy;
    mp->xlower = xlower;
    mp->ylower = ylower;
    mp->xupper = xupper;
    mp->yupper = yupper;

    /* Set up area for storage - this is needed for ghost patches, 
    and updated patches */
    {
        int ll[2];
        int ur[2];
        for (int idir = 0; idir < 2; idir++)
        {
            ll[idir] = -mbc;
        }
        ur[0] = mx + mbc + 1;
        ur[1] = my + mbc + 1;

        Box box_p(ll,ur);
        mp->area.define(box_p,1);
    }


    if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED
        && build_mode == FCLAW2D_BUILD_FOR_UPDATE)
    {
        /* 
        From here on out, we build for patches that get updated 

        24 additional field variables needed for all metric terms
            xp,yp,zp           : 3
            xd,yd,zd           : 3
            surf_normals       : 3
            curvature          : 1
            <xy>face normals   : 6
            <xy>face tangents  : 6
            edge lengths       : 2
            -----------------------
            Total              : 24

            We should come up with a way to store only what is needed 
        */

        int ll[2];
        int ur[2];

        for (int idir = 0; idir < 2; idir++)
        {
            ll[idir] = -mbc;
        }
        ur[0] = mx + mbc + 1;
        ur[1] = my + mbc + 1;

        Box box_p(ll,ur);   /* Store cell centered values here */

        /* Mesh cell centers of physical mesh */
        mp->xp.define(box_p,1);
        mp->yp.define(box_p,1);
        mp->zp.define(box_p,1);
        mp->surf_normals.define(box_p,3);
        mp->curvature.define(box_p,1);

        /* Node centered values */
        for (int idir = 0; idir < 2; idir++)
        {
            ll[idir] = -mbc;
        }
        ur[0] = mx + mbc + 2;
        ur[1] = my + mbc + 2;
        Box box_d(ll,ur);

        mp->xd.define(box_d,1);
        mp->yd.define(box_d,1);
        mp->zd.define(box_d,1);

        /* Face centered values */
        mp->xface_normals.define(box_d,3);
        mp->yface_normals.define(box_d,3);
        mp->xface_tangents.define(box_d,3);
        mp->yface_tangents.define(box_d,3);
        mp->edge_lengths.define(box_d,2);
    }
}


static
void metric_average_area_from_fine(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *fine_patches,
                                   fclaw2d_patch_t *coarse_patch,
                                   int blockno, 
                                   int coarse_patchno,
                                   int fine0_patchno)

{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt();
    int mx,my, mbc;
    double xlower,ylower,dx,dy;

    double *areacoarse, *areafine;
    int igrid;

    fclaw2d_metric_patch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                   &xlower,&ylower,&dx,&dy);

    areacoarse = fclaw2d_metric_patch_get_area(coarse_patch);

    for(igrid = 0; igrid < 4; igrid++)
    {
        areafine = fclaw2d_metric_patch_get_area(&fine_patches[igrid]);

        FCLAW2D_FORT_AVERAGE_AREA(&mx,&my,&mbc,areacoarse,areafine,&igrid);
    }

    metric_vt->compute_area_ghost(glob,coarse_patch,blockno,coarse_patchno);
}

/* --------------------------------- Public interface  -------------------------------- */

void fclaw2d_metric_patch_compute_area (fclaw2d_global_t *glob,
                                       fclaw2d_patch_t* this_patch,
                                       int blockno, int patchno)
{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt();
    FCLAW_ASSERT(metric_vt->compute_area);

    metric_vt->compute_area(glob,this_patch,blockno,patchno);
}


void fclaw2d_metric_patch_setup(fclaw2d_global_t* glob,
                                fclaw2d_patch_t* this_patch,
                                int blockno,
                                int patchno)
{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt();

    /* Setup the mesh and normals */
    metric_vt->compute_mesh   (glob,this_patch,blockno,patchno);
    metric_vt->compute_normals(glob,this_patch,blockno,patchno);
}



void fclaw2d_metric_patch_setup_from_fine(fclaw2d_global_t *glob,
                                          fclaw2d_patch_t *fine_patches,
                                          fclaw2d_patch_t *coarse_patch,
                                          int blockno,
                                          int coarse_patchno,
                                          int fine0_patchno)
{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt();

    metric_average_area_from_fine(glob,fine_patches,coarse_patch,
                                  blockno, coarse_patchno, 
                                  fine0_patchno);

    metric_vt->compute_mesh(glob,coarse_patch,blockno,coarse_patchno);
    metric_vt->compute_normals(glob,coarse_patch,blockno,coarse_patchno);
}



/* ------------------------------------ Virtual table  -------------------------------- */

static
fclaw2d_metric_vtable_t* metric_vt_init()
{
    FCLAW_ASSERT(s_metric_vt.is_set == 0);
    return &s_metric_vt;
}

fclaw2d_metric_vtable_t* fclaw2d_metric_vt()
{
    FCLAW_ASSERT(s_metric_vt.is_set != 0);
    return &s_metric_vt;
}


void fclaw2d_metric_vtable_initialize()  
{
    fclaw2d_metric_vtable_t *metric_vt = metric_vt_init();

    metric_vt->compute_mesh          = fclaw2d_metric_compute_mesh_default;
    metric_vt->compute_area          = fclaw2d_metric_compute_area_default;
    metric_vt->compute_area_ghost    = fclaw2d_metric_compute_area_ghost_default;
    metric_vt->compute_normals       = fclaw2d_metric_compute_normals_default;

    /* Fortran files */
    metric_vt->fort_compute_mesh          = &FCLAW2D_FORT_COMPUTE_MESH;
    metric_vt->fort_compute_normals       = &FCLAW2D_FORT_COMPUTE_NORMALS;
    metric_vt->fort_compute_tangents      = &FCLAW2D_FORT_COMPUTE_TANGENTS;
    metric_vt->fort_compute_surf_normals  = &FCLAW2D_FORT_COMPUTE_SURF_NORMALS;

    metric_vt->is_set = 1;
}


/* --------------------------------- Misc access functions ---------------------------- */

/* These functions are not virtualized and are not defined by the 
   patch interface */

fclaw2d_metric_patch_t* fclaw2d_metric_get_metric_patch(fclaw2d_patch_t* this_patch)

{
    return get_metric_patch(this_patch);
}

double* fclaw2d_metric_patch_get_area(fclaw2d_patch_t* this_patch)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);
    return mp->area.dataPtr();
}



void fclaw2d_metric_patch_scalar(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch,
                                 double **area, double** edgelengths,
                                 double **curvature)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);
    *area = mp->area.dataPtr();
    *edgelengths =  mp->edge_lengths.dataPtr();
    *curvature = mp->curvature.dataPtr();
}


void fclaw2d_metric_patch_vector(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch,
                                 double **xnormals, double **ynormals,
                                 double **xtangents, double **ytangents,
                                 double **surfnormals)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);
    *xnormals = mp->xface_normals.dataPtr();
    *ynormals = mp->yface_normals.dataPtr();
    *xtangents = mp->xface_tangents.dataPtr();
    *ytangents = mp->yface_tangents.dataPtr();
    *surfnormals = mp->surf_normals.dataPtr();
}


void fclaw2d_metric_patch_grid_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    int* mx, int* my, int* mbc,
                                    double* xlower, double* ylower,
                                    double* dx, double* dy)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);
    *mx     = mp->mx;
    *my     = mp->my;
    *mbc    = mp->mbc;
    *xlower = mp->xlower;
    *ylower = mp->ylower;
    *dx     = mp->dx;
    *dy     = mp->dy;
}


void fclaw2d_metric_patch_mesh_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    double **xp, double **yp, double **zp,
                                    double **xd, double **yd, double **zd,
                                    double **area)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);
    *xp = mp->xp.dataPtr();
    *yp = mp->yp.dataPtr();
    *zp = mp->zp.dataPtr();
    *xd = mp->xd.dataPtr();
    *yd = mp->yd.dataPtr();
    *zd = mp->zd.dataPtr();
    *area = mp->area.dataPtr();
}

void fclaw2d_metric_patch_mesh_data2(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch,
                                     double **xnormals, double **ynormals,
                                     double **xtangents, double **ytangents,
                                     double **surfnormals,
                                     double **edgelengths, double **curvature)
{
    fclaw2d_metric_patch_t* mp = get_metric_patch(this_patch);
    *xnormals    = mp->xface_normals.dataPtr();
    *ynormals    = mp->yface_normals.dataPtr();
    *xtangents   = mp->xface_tangents.dataPtr();
    *ytangents   = mp->yface_tangents.dataPtr();
    *surfnormals = mp->surf_normals.dataPtr();
    *curvature   = mp->curvature.dataPtr();
    *edgelengths = mp->edge_lengths.dataPtr();
}






