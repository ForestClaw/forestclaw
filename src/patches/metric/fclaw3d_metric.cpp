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


#if 0
/* Fix syntax highlighting */
#endif

#define METRIC_VTABLE_NAME "fclaw3d_metric"

#include "fclaw3d_metric.h"
#include "fclaw3d_metric.hpp"

#include <fclaw_pointer_map.h>

#include <fclaw_global.h>
#include <fclaw_patch.h>  
#include <fclaw3d_defs.h>

static
fclaw3d_metric_patch_t* get_metric_patch(fclaw_global_t* glob,
                                         fclaw_patch_t *patch)
{
    return (fclaw3d_metric_patch_t*) fclaw_patch_metric_patch(glob, patch);
}


/* ----------------------------- Creating/deleting patches ---------------------------- */

fclaw3d_metric_patch_t* fclaw3d_metric_patch_new()
{
    fclaw3d_metric_patch_t *mp = new fclaw3d_metric_patch_t;
    return mp;
}

void fclaw3d_metric_patch_delete(fclaw3d_metric_patch_t **mp)
{
    FCLAW_ASSERT(mp != NULL);
    delete *mp;
    *mp = NULL;
}


void fclaw3d_metric_patch_define(fclaw_global_t* glob,
                                 fclaw_patch_t* patch, 
                                 int mx, int my, int mz, int mbc, 
                                 double dx, double dy,double dz,
                                 double xlower, double ylower, double zlower,
                                 double xupper, double yupper, double zupper,
                                 int blockno, int patchno, 
                                 fclaw_build_mode_t build_mode)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);

    mp->mx   = mx;
    mp->my   = my;
    mp->mz   = mz;
    mp->mbc  = mbc;
    mp->blockno = blockno;

    mp->dx = dx;
    mp->dy = dy;
    mp->dz = dz;
    mp->xlower = xlower;
    mp->ylower = ylower;
    mp->zlower = zlower;
    mp->xupper = xupper;
    mp->yupper = yupper;
    mp->zupper = zupper;


    /* Set up area for storage - this is needed for ghost patches, 
    and updated patches */
    Box box_p, box_d;
    {
        /* Create box for primary grid (cell-centered) */
        int ll_p[FCLAW3D_SPACEDIM];
        int ur_p[FCLAW3D_SPACEDIM];
        for (int idir = 0; idir < FCLAW3D_SPACEDIM; idir++)
        {
            ll_p[idir] = -mbc;
        }
        ur_p[0] = mp->mx + mbc + 1;
        ur_p[1] = mp->my + mbc + 1;
        ur_p[2] = mp->mz + mbc + 1;

        box_p = Box(ll_p,ur_p,FCLAW3D_SPACEDIM);


        /* Create box for nodes (cell corners) */
        int ll_d[FCLAW3D_SPACEDIM];
        int ur_d[FCLAW3D_SPACEDIM];

        /* Create node centered box */
        for (int idir = 0; idir < FCLAW3D_SPACEDIM; idir++)
        {
            ll_d[idir] = -mbc;
        }
        ur_d[0] = mp->mx + mbc + 2;
        ur_d[1] = mp->my + mbc + 2;
        ur_d[2] = mp->mz + mbc + 2;

        box_d = Box(ll_d,ur_d,FCLAW3D_SPACEDIM);


        /* volume is needed in ghost patches for averaging.   The face area is 
           only needed if we are updating a cell, but since volume and 
           face area are computed together, we should also allocate memory 
           for face areas
        */
        mp->volume.define(box_p,1);
        mp->face_area.define(box_d,3);

        /* Node values are needed by parallel ghost patches when we have an 
           affine map.  In this case, the affine map can use pre-computed xd,yd,zd
           values.

           Note: This should be revisited since we are doing extra work for the affine
           map. 
        */
        mp->xd.define(box_d,1);
        mp->yd.define(box_d,1);
        mp->zd.define(box_d,1);
    }


    /* Only allocate memory that is needed */
    if (build_mode != FCLAW_BUILD_FOR_GHOST_AREA_PACKED
       && build_mode == FCLAW_BUILD_FOR_UPDATE)    
    {
        /* 
        From here on out, we build for patches that get updated 

        In 2d manifold : 
        24 additional field variables needed for all metric terms on a 
        manifold
            xp,yp,zp           : 3
            xd,yd,zd           : 3
            surf_normals       : 3    (3x1 vector)
            curvature          : 1
            area               : 1
            <xy>face normals   : 6    (2 3x1 vectors)
            <xy>face tangents  : 6    (2 3x1 vectors)
            edge lengths       : 2
            -----------------------
            Total              : 25

        In 3d : 
        36 additional field variables needed for all metric terms
            xp,yp,zp           : 3
            xd,yd,zd           : 3
            face areas         : 3
            volume             : 1
            rotation matrix    : 27 (3 3x3 matrices)
            -----------------------
            Total              : 37

            We should come up with a way to store only what is needed 
        */

        /* Mesh cell centers of physical mesh */
        mp->xp.define(box_p,1);
        mp->yp.define(box_p,1);
        mp->zp.define(box_p,1);

        /* Store face areas of left, front, bottom edge of box */
        //mp->face_area.define(box_d,3);

        /* Store 3x3 rotation matrix for each of three faces */
        mp->xrot.define(box_d,9);
        mp->yrot.define(box_d,9);
        mp->zrot.define(box_d,9);

    }
}



/* For 3d extruded mesh, this averages cell volumes */
static
void metric_average_area_from_fine(fclaw_global_t *glob,
                                   fclaw_patch_t *fine_patches,
                                   fclaw_patch_t *coarse_patch,
                                   int blockno, 
                                   int coarse_patchno,
                                   int fine0_patchno)
{
    if(glob->domain->refine_dim == 3)
    {
        fclaw_abortf("Metric averaging not implemented for oct-trees.\n");
    }

    fclaw3d_metric_vtable_t *metric_vt = fclaw3d_metric_vt(glob);
    int mx,my,mz, mbc;
    double xlower,ylower,zlower,dx,dy,dz;

    fclaw3d_metric_patch_grid_data(glob,coarse_patch,&mx,&my,&mz,&mbc,
                                   &xlower,&ylower,&zlower,&dx,&dy,&dz);    
    double *volcoarse, *fa_coarse;
    fclaw3d_metric_patch_scalar(glob,coarse_patch,&volcoarse, &fa_coarse);


    for(int igrid = 0; igrid < 4; igrid++)
    {
        double *volfine, *fa_fine;
        fclaw3d_metric_patch_scalar(glob,&fine_patches[igrid],&volfine, &fa_fine);

        /* Average from fine to coarse when creating coarse grid from a fine grid. 
           This will be exact, since both volume and face areas are computing at 
           finest level resolution.
        */
        FCLAW3DX_METRIC_FORT_AVERAGE_VOLUME(&mx,&my,&mz, &mbc,
                                            volcoarse,volfine, 
                                            &igrid);

        FCLAW3DX_METRIC_FORT_AVERAGE_FACEAREA(&mx,&my,&mz, &mbc,
                                              fa_coarse,fa_fine, &igrid);
    }

/* Compute volume from scratch in second layer of coarse grid ghost cells. These
   are not averaged from the fine grid (first layer can be averaged from fine
   grid, however)
*/
    metric_vt->compute_volume_ghost(glob,coarse_patch,blockno,coarse_patchno);
}

/* --------------------------------- Public interface  -------------------------------- */

void fclaw3d_metric_patch_compute_volume (fclaw_global_t *glob,
                                          fclaw_patch_t* patch,
                                          int blockno, int patchno)
{
    fclaw3d_metric_vtable_t *metric_vt = fclaw3d_metric_vt(glob);
    FCLAW_ASSERT(metric_vt->compute_volume);
    metric_vt->compute_volume(glob,patch,blockno,patchno);
}


void fclaw3d_metric_patch_build(fclaw_global_t* glob,
                                fclaw_patch_t* patch,
                                int blockno,
                                int patchno)
{
    fclaw3d_metric_vtable_t *metric_vt = fclaw3d_metric_vt(glob);
    FCLAW_ASSERT(metric_vt != NULL);

    /* Compute (xp,yp,zp) and (xd,yd,zd) */
    metric_vt->compute_mesh(glob,patch,blockno,patchno);

    /* Compute areas/volumes ($$$) from scratch. 
       Note : These are all computed on finest level 
       mesh and averaged down to coarser meshes.  This is 
       required from geometric consistency */
    /* Compute 3d volumes and 2d face areas */
    metric_vt->compute_volume(glob,patch,blockno,patchno);
    if (metric_vt->compute_basis != NULL)
    {
        /* In 2d : Surface normals, tangents, edge lengths, 
           surface normals and curvature. 

           In 3d : Rotation matrix at each face. 
        */
        metric_vt->compute_basis(glob,patch,blockno,patchno);
    }
}



void fclaw3d_metric_patch_build_from_fine(fclaw_global_t *glob,
                                          fclaw_patch_t *fine_patches,
                                          fclaw_patch_t *coarse_patch,
                                          int blockno,
                                          int coarse_patchno,
                                          int fine0_patchno)
{
    /* This routine does a complete build using fine grid data to average
       volumes and in 3d, face areas */
    fclaw3d_metric_vtable_t *metric_vt = fclaw3d_metric_vt(glob);

    /* Compute xd,yd,zd, xp,yp,zp */
    metric_vt->compute_mesh(glob,coarse_patch,blockno,coarse_patchno);

    /* Compute areas/volumes by averaging from finer grids  */

    metric_average_area_from_fine(glob,fine_patches,coarse_patch,
                                  blockno, coarse_patchno, 
                                  fine0_patchno);

    if (metric_vt->compute_basis != NULL)
    {
        /* In 2d : Surface normals and tangents at each face
           In 3d : Rotation matrix for each face. 

           Note : These are not averaged from finer grids, but are 
           built from scratch here. 
        */
        metric_vt->compute_basis(glob,coarse_patch,blockno,coarse_patchno);        
    }
}



/* ------------------------------------ Virtual table  -------------------------------- */

static
fclaw3d_metric_vtable_t* metric_vt_new()
{
    return (fclaw3d_metric_vtable_t*) FCLAW_ALLOC_ZERO (fclaw3d_metric_vtable_t, 1);
}

static
void metric_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

fclaw3d_metric_vtable_t* fclaw3d_metric_vt(fclaw_global_t* glob)
{
	fclaw3d_metric_vtable_t* metric_vt = (fclaw3d_metric_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, METRIC_VTABLE_NAME);
	FCLAW_ASSERT(metric_vt != NULL);
	FCLAW_ASSERT(metric_vt->is_set != 0);
	return metric_vt;
}


void fclaw3d_metric_vtable_initialize(fclaw_global_t* glob)  
{
    fclaw3d_metric_vtable_t *metric_vt = metric_vt_new();


    /* Fortran files */
    metric_vt->compute_mesh          = fclaw3d_metric_compute_mesh_default;
    metric_vt->compute_volume        = fclaw3d_metric_compute_volume_default;
    metric_vt->compute_volume_ghost  = fclaw3d_metric_compute_volume_ghost_default;
    metric_vt->compute_basis         = fclaw3d_metric_compute_basis_default;    

    metric_vt->fort_compute_mesh     = &FCLAW3D_METRIC_FORT_COMPUTE_MESH;
    metric_vt->fort_compute_volume   = &FCLAW3D_METRIC_FORT_COMPUTE_VOLUME;
    metric_vt->fort_compute_basis    = &FCLAW3D_METRIC_FORT_COMPUTE_BASIS;

    metric_vt->is_set = 1;

	if(fclaw_pointer_map_get(glob->vtables,METRIC_VTABLE_NAME) != NULL)
    {
        fclaw_abortf("Metric vtable %s already set\n",METRIC_VTABLE_NAME);
    }
	fclaw_pointer_map_insert(glob->vtables,METRIC_VTABLE_NAME, metric_vt, metric_vt_destroy);
}


/* --------------------------------- Misc access functions ---------------------------- */

/* These functions are not virtualized and are not defined by the 
   patch interface */

fclaw3d_metric_patch_t* fclaw3d_metric_get_metric_patch(fclaw_global_t* glob,
                                                        fclaw_patch_t* patch)

{
    return get_metric_patch(glob, patch);
}

double* fclaw3d_metric_patch_get_volume(fclaw_global_t* glob,
                                        fclaw_patch_t* patch)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);
    return mp->volume.dataPtr();
}


/* ----------- See fclaw3d_metric.cpp for 3d versions of functions below -------------- */

void fclaw3d_metric_patch_scalar(fclaw_global_t* glob,
                                 fclaw_patch_t* patch,
                                 double **volume, double** faceareas)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);
    *volume = mp->volume.dataPtr();
    *faceareas =  mp->face_area.dataPtr();
}



void fclaw3d_metric_patch_basis(fclaw_global_t* glob,
                                fclaw_patch_t* patch,
                                double **xrot, double **yrot, double **zrot)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);
    *xrot = mp->xrot.dataPtr();
    *yrot = mp->yrot.dataPtr();
    *zrot = mp->zrot.dataPtr();
}

void fclaw3d_metric_patch_grid_data(fclaw_global_t* glob,
                                    fclaw_patch_t* patch,
                                    int* mx, int* my, int* mz, 
                                    int* mbc,
                                    double* xlower, double* ylower, double *zlower,
                                    double* dx, double* dy, double* dz)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);
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


void fclaw3d_metric_patch_mesh_data(fclaw_global_t* glob,
                                    fclaw_patch_t* patch,
                                    double **xp, double **yp, double **zp,
                                    double **xd, double **yd, double **zd,
                                    double **volume, double** faceareas)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);
    *xp = mp->xp.dataPtr();
    *yp = mp->yp.dataPtr();
    *zp = mp->zp.dataPtr();
    *xd = mp->xd.dataPtr();
    *yd = mp->yd.dataPtr();
    *zd = mp->zd.dataPtr();
    *volume = mp->volume.dataPtr();
    *faceareas = mp->face_area.dataPtr();
}

int fclaw3d_metric_patch_nodes_size(fclaw_global_t* glob,
                                    fclaw_patch_t* patch)
{
    fclaw3d_metric_patch_t* mp = get_metric_patch(glob, patch);
    return mp->xd.size();
}