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
#include <fclaw_patch.h>

#include <fclaw_clawpatch_pillow.h>

#include <fclaw_clawpatch.h>
#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>
#include <fclaw3dx_clawpatch46_fort.h>


#include <fclaw3d_metric.h>

#include <fclaw_pointer_map.h>

struct fclaw_patch_transform_data;  /* Not used here, so we leave it incomplete */

static
void pillow_copy_block_corner(fclaw_global_t* glob,
                              fclaw_patch_t* patch, 
                              fclaw_patch_t *corner_patch,
                              int blockno,
                              int corner_blockno,
                              int icorner,
                              int time_interp,
                              struct fclaw_patch_transform_data *transform_data)
{
    int meqn;
    double *qthis;    
    fclaw_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);

    double *qcorner = fclaw_clawpatch_get_q(glob,corner_patch);

    fclaw_clawpatch_pillow_vtable_t* pillow_vt = fclaw_clawpatch_pillow_vt(glob);
    FCLAW_ASSERT(pillow_vt != NULL);

    if(pillow_vt->dim == 2)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;
        fclaw_clawpatch_grid_data_2d(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        pillow_vt->d2->fort_copy_block_corner(&mx, &my, &mbc, &meqn, 
                                              qthis, qcorner,
                                              &icorner, &blockno);
    }
    else 
    {
        int mx,my,mz,mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_grid_data_3d(glob,patch,&mx,&my,&mz,&mbc,
                                    &xlower,&ylower,&zlower, &dx, &dy, &dz);
        FCLAW_ASSERT(pillow_vt->d3->fort_copy_block_corner != NULL);

        pillow_vt->d3->fort_copy_block_corner(&mx, &my, &mz, &mbc, &meqn, 
                                              qthis, qcorner,
                                              &icorner, &blockno);
    }

}

static
void pillow_average_block_corner(fclaw_global_t *glob,
                                 fclaw_patch_t* coarse_patch,
                                 fclaw_patch_t *fine_patch,
                                 int coarse_blockno,
                                 int fine_blockno,
                                 int icorner_coarse,
                                 int time_interp,
                                 struct fclaw_patch_transform_data* transform_data)
{

    int refratio = 2;

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,
                                    &qcoarse,&meqn);
    double* qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    double *areacoarse = fclaw2d_clawpatch_get_area(glob,coarse_patch);
    double *areafine = fclaw2d_clawpatch_get_area(glob,fine_patch);

    fclaw_clawpatch_pillow_vtable_t* pillow_vt = fclaw_clawpatch_pillow_vt(glob);

    if(pillow_vt->dim == 2)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;

        fclaw_clawpatch_grid_data_2d(glob,coarse_patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);


        pillow_vt->d2->fort_average_block_corner(&mx,&my,&mbc,&meqn,
                                                 &refratio,qcoarse,qfine,
                                                 areacoarse,areafine,
                                                 &icorner_coarse,&coarse_blockno);
    }
    else
    {
        int mx,my,mz,mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_grid_data_3d(glob,coarse_patch,&mx,&my,&mz, &mbc,
                                    &xlower,&ylower, &zlower, &dx,&dy, &dz);

        pillow_vt->d3->fort_average_block_corner(&mx,&my,&mz, &dz, &mbc,&meqn,
                                                 &refratio,qcoarse,qfine,
                                                 areacoarse,areafine,
                                                 &icorner_coarse,&coarse_blockno);
    }
}

static
void pillow_interpolate_block_corner(fclaw_global_t* glob,
                                     fclaw_patch_t* coarse_patch,
                                     fclaw_patch_t *fine_patch,
                                     int coarse_blockno,
                                     int fine_blockno,
                                     int icoarse_corner,
                                     int time_interp,
                                     struct fclaw_patch_transform_data* transform_data)

{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,
                                    &qcoarse,&meqn);

    double* qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    fclaw_clawpatch_pillow_vtable_t* pillow_vt = fclaw_clawpatch_pillow_vt(glob);
    int refratio = 2;
    if(pillow_vt->dim == 2)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;
        fclaw_clawpatch_grid_data_2d(glob,coarse_patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        pillow_vt->d2->fort_interpolate_block_corner(&mx, &my, &mbc, &meqn,
                                                 &refratio, qcoarse, qfine,
                                                 &icoarse_corner, &coarse_blockno);
    }
    else 
    {
        int mx,my,mz,mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_grid_data_3d(glob,coarse_patch,&mx,&my,&mz, &mbc,
                                    &xlower,&ylower,&zlower, &dx,&dy, &dz);

        pillow_vt->d3->fort_interpolate_block_corner(&mx, &my, &mz, &mbc, &meqn,
                                                     &refratio, qcoarse, qfine,
                                                     &icoarse_corner, &coarse_blockno);
    }
}

/* ----------------------------- Use pillow sphere ------------------------------------ */

void fclaw_clawpatch_use_pillowsphere(fclaw_global_t* glob)
{
    fclaw_patch_vtable_t* patch_vt = fclaw_patch_vt(glob);

    patch_vt->d2->copy_block_corner          = pillow_copy_block_corner;
    patch_vt->d2->average_block_corner       = pillow_average_block_corner;
    patch_vt->d2->interpolate_block_corner   = pillow_interpolate_block_corner;

}


/* -------------------------------- Virtual table ------------------------------------- */

static
fclaw_clawpatch_pillow_vtable_t* pillow_vt_new(int dim)
{
    fclaw_clawpatch_pillow_vtable_t* pillow_vt =
           FCLAW_ALLOC_ZERO (fclaw_clawpatch_pillow_vtable_t, 1);
    pillow_vt->dim = dim;
    if(dim == 2)
    {
        pillow_vt->d2 = FCLAW_ALLOC_ZERO(struct fclaw2d_clawpatch_pillow_vtable,1);
    }
    else 
    {
        pillow_vt->d3 = FCLAW_ALLOC_ZERO(struct fclaw3d_clawpatch_pillow_vtable,1);
    }
    return pillow_vt;
}

static
void pillow_vt_destroy(void* vt)
{
    fclaw_clawpatch_pillow_vtable_t *pillow_vt = (fclaw_clawpatch_pillow_vtable_t*) vt;
    if (pillow_vt->dim == 2)
    {
        FCLAW_FREE (pillow_vt->d2);
    }
    else
    {
        FCLAW_FREE (pillow_vt->d3);
    }
    FCLAW_FREE (pillow_vt);
}

static
void fclaw_clawpatch_pillow_vtable_initialize(int dim, 
                                              fclaw_global_t* glob,
                                              int claw_version)
{
    fclaw_clawpatch_pillow_vtable_t *pillow_vt = pillow_vt_new(dim);

    if(dim == 2)
    {
        if (claw_version == 4)
        {
            pillow_vt->d2->fort_copy_block_corner        = FCLAW2D_PILLOW46_COPY_BLOCK_CORNER;
            pillow_vt->d2->fort_average_block_corner     = FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER;
            pillow_vt->d2->fort_interpolate_block_corner = FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER;
        }
        else if (claw_version == 5)
        {
            pillow_vt->d2->fort_copy_block_corner        = FCLAW2D_PILLOW5_COPY_BLOCK_CORNER;
            pillow_vt->d2->fort_average_block_corner     = FCLAW2D_PILLOW5_AVERAGE_BLOCK_CORNER;
            pillow_vt->d2->fort_interpolate_block_corner = FCLAW2D_PILLOW5_INTERPOLATE_BLOCK_CORNER;
        }
    }
    else
    {
        if (claw_version == 4)
        {
            pillow_vt->d3->fort_copy_block_corner        = FCLAW3DX_PILLOW46_COPY_BLOCK_CORNER;
            pillow_vt->d3->fort_average_block_corner     = FCLAW3DX_PILLOW46_AVERAGE_BLOCK_CORNER;
            pillow_vt->d3->fort_interpolate_block_corner = FCLAW3DX_PILLOW46_INTERPOLATE_BLOCK_CORNER;
        }
    }

    pillow_vt->is_set = 1;

    FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables, "fclaw_clawpatch_pillow_vtable") == NULL);
    fclaw_pointer_map_insert(glob->vtables, "fclaw_clawpatch_pillow_vtable", pillow_vt, pillow_vt_destroy);

}

void fclaw2d_clawpatch_pillow_vtable_initialize(fclaw_global_t* glob,
                                                int claw_version)
{
    fclaw_clawpatch_pillow_vtable_initialize(2,glob,claw_version);
}

void fclaw3d_clawpatch_pillow_vtable_initialize(fclaw_global_t* glob,
                                                int claw_version)
{
    fclaw_clawpatch_pillow_vtable_initialize(3,glob,claw_version);
}


/* ------------------------------- Public access functions ---------------------------- */


fclaw_clawpatch_pillow_vtable_t* fclaw_clawpatch_pillow_vt(fclaw_global_t* glob)
{

    fclaw_clawpatch_pillow_vtable_t* pillow_vt = (fclaw_clawpatch_pillow_vtable_t*) 
        fclaw_pointer_map_get(glob->vtables, 
                              "fclaw_clawpatch_pillow_vtable");

    FCLAW_ASSERT(pillow_vt != NULL);
    FCLAW_ASSERT(pillow_vt->is_set != 0);
    return pillow_vt;
}