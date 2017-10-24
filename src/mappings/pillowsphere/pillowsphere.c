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


#include <pillowsphere.h>

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_transform.h>

void pillow_copy_block_corner(fclaw2d_global_t* glob,
                              fclaw2d_patch_t* this_patch, 
                              fclaw2d_patch_t *corner_patch,
                              int icorner,
                              int time_interp,
                              fclaw2d_patch_transform_data_t *transform_data)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qthis,*qcorner;

    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);

    qcorner = fclaw2d_clawpatch_get_q(glob,corner_patch);

    int blockno = fclaw2d_patch_get_blockno(this_patch);
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    FCLAW2D_FORT_MB_EXCHANGE_BLOCK_CORNER_GHOST(&mx, &my, &mbc, &meqn, 
                                                qthis, qcorner,
                                                &icorner, &blockno);

}

void pillow_average_block_corner(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t* coarse_patch,
                                 fclaw2d_patch_t *fine_patch,
                                 int icorner_coarse,
                                 int time_interp,
                                 fclaw2d_patch_transform_data_t* transform_data)
{
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

    int blockno = fclaw2d_patch_get_blockno(coarse_patch);
    FCLAW2D_FORT_MB_AVERAGE_BLOCK_CORNER_GHOST(&mx,&my,&mbc,&meqn,
                                               &refratio,qcoarse,qfine,
                                               areacoarse,areafine,
                                               &icorner_coarse,&blockno);
}

void pillow_interpolate_block_corner(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* coarse_patch,
                                     fclaw2d_patch_t *fine_patch,
                                     int icoarse_corner,
                                     int time_interp,
                                     fclaw2d_patch_transform_data_t* transform_data)

{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qcoarse, *qfine;

    int refratio = 2;

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,
                                    &qcoarse,&meqn);

    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    fclaw2d_clawpatch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int blockno = fclaw2d_patch_get_blockno(coarse_patch);

    FCLAW2D_FORT_MB_INTERPOLATE_BLOCK_CORNER_GHOST(&mx, &my, &mbc, &meqn,
                                                   &refratio, qcoarse, qfine,
                                                   &icoarse_corner, &blockno);
}




