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


#include <fclaw2d_clawpatch_pillow.h>

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>

static fclaw2d_clawpatch_pillow_vtable_t s_clawpatch_pillow_vt;

struct fclaw2d_patch_transform_data;  /* Not used here, so we leave it incomplete */

static
void pillow_copy_block_corner(fclaw2d_global_t* glob,
                              fclaw2d_patch_t* this_patch, 
                              fclaw2d_patch_t *corner_patch,
                              int this_blockno,
                              int corner_blockno,
                              int icorner,
                              int time_interp,
                              struct fclaw2d_patch_transform_data *transform_data)
{
    fclaw2d_clawpatch_pillow_vtable_t* pillow_vt = fclaw2d_clawpatch_pillow_vt();

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qthis,*qcorner;

    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);

    qcorner = fclaw2d_clawpatch_get_q(glob,corner_patch);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

#if FCLAW2D_PATCHDIM == 2
    pillow_vt->fort_copy_block_corner(&mx, &my, &mbc, &meqn, 
                                      qthis, qcorner,
                                      &icorner, &this_blockno);
#endif

}

static
void pillow_average_block_corner(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t* coarse_patch,
                                 fclaw2d_patch_t *fine_patch,
                                 int coarse_blockno,
                                 int fine_blockno,
                                 int icorner_coarse,
                                 int time_interp,
                                 struct fclaw2d_patch_transform_data* transform_data)
{
    fclaw2d_clawpatch_pillow_vtable_t* pillow_vt = fclaw2d_clawpatch_pillow_vt();

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qcoarse, *qfine;

    int refratio = 2;

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,
                                    &qcoarse,&meqn);

    double *areacoarse = fclaw2d_clawpatch_get_area(glob,coarse_patch);
    double *areafine = fclaw2d_clawpatch_get_area(glob,fine_patch);

    fclaw2d_clawpatch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

#if FCLAW2D_PATCHDIM == 2
    pillow_vt->fort_average_block_corner(&mx,&my,&mbc,&meqn,
                                         &refratio,qcoarse,qfine,
                                         areacoarse,areafine,
                                         &icorner_coarse,&coarse_blockno);
#endif
}

static
void pillow_interpolate_block_corner(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* coarse_patch,
                                     fclaw2d_patch_t *fine_patch,
                                     int coarse_blockno,
                                     int fine_blockno,
                                     int icoarse_corner,
                                     int time_interp,
                                     struct fclaw2d_patch_transform_data* transform_data)

{
    fclaw2d_clawpatch_pillow_vtable_t* pillow_vt = fclaw2d_clawpatch_pillow_vt();

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qcoarse, *qfine;

    int refratio = 2;

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,
                                    &qcoarse,&meqn);

    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    fclaw2d_clawpatch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

#if FCLAW2D_PATCHDIM == 2
    pillow_vt->fort_interpolate_block_corner(&mx, &my, &mbc, &meqn,
                                             &refratio, qcoarse, qfine,
                                             &icoarse_corner, &coarse_blockno);
#endif
}

/* ----------------------------- Use pillow sphere ------------------------------------ */

void fclaw2d_clawpatch_use_pillowsphere()
{
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();

    patch_vt->copy_block_corner          = pillow_copy_block_corner;
    patch_vt->average_block_corner       = pillow_average_block_corner;
    patch_vt->interpolate_block_corner   = pillow_interpolate_block_corner;

}


/* -------------------------------- Virtual table ------------------------------------- */

static
fclaw2d_clawpatch_pillow_vtable_t* pillow_vt_init()
{
    FCLAW_ASSERT(s_clawpatch_pillow_vt.is_set == 0);
    return &s_clawpatch_pillow_vt;
}

void fclaw2d_clawpatch_pillow_vtable_initialize(int claw_version)
{
    fclaw2d_clawpatch_pillow_vtable_t *pillow_vt = pillow_vt_init();

#if FCLAW2D_PATCHDIM == 2
    if (claw_version == 4)
    {
        pillow_vt->fort_copy_block_corner        = FCLAW2D_PILLOW46_COPY_BLOCK_CORNER;
        pillow_vt->fort_average_block_corner     = FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER;
        pillow_vt->fort_interpolate_block_corner = FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER;
    }
    else if (claw_version == 5)
    {
        pillow_vt->fort_copy_block_corner        = FCLAW2D_PILLOW5_COPY_BLOCK_CORNER;
        pillow_vt->fort_average_block_corner     = FCLAW2D_PILLOW5_AVERAGE_BLOCK_CORNER;
        pillow_vt->fort_interpolate_block_corner = FCLAW2D_PILLOW5_INTERPOLATE_BLOCK_CORNER;
    }
#endif

    pillow_vt->is_set = 1;
}


/* ------------------------------- Public access functions ---------------------------- */

fclaw2d_clawpatch_pillow_vtable_t* fclaw2d_clawpatch_pillow_vt()
{
    FCLAW_ASSERT(s_clawpatch_pillow_vt.is_set != 0);
    return &s_clawpatch_pillow_vt;
}


