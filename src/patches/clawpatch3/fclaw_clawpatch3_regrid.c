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

#include <fclaw2d_forestclaw.h>
#include <fclaw_clawpatch3.h>
#include <fclaw2d_metric_default_fort.h>

int fclaw_clawpatch3_tag4refinement(fclaw2d_global_t *glob,
                                     fclaw2d_patch_t *this_patch,
                                     int blockno, int patchno,
                                     int initflag)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);

    int mx,my,mz,mbc,meqn;
    double xlower,ylower,zlower,dx,dy,dz;
    double *q;
    int tag_patch;
    double refine_threshold;

    refine_threshold = gparms->refine_threshold;

    fclaw_clawpatch3_grid_data(glob,this_patch,&mx,&my,&mz,&mbc,
                                 &xlower,&ylower,&zlower,&dx,&dy,&dz);

    fclaw_clawpatch3_soln_data(glob,this_patch,&q,&meqn);

    tag_patch = 0;
    fclaw_clawpatch3_vt()->fort_tag4refinement(&mx,&my,&mz,&mbc,&meqn,&xlower,&ylower,&zlower,
                                                 &dx,&dy,&dz,
                                                 &blockno, q,&refine_threshold,
                                                 &initflag,&tag_patch);
    return tag_patch;
}

int fclaw_clawpatch3_tag4coarsening(fclaw2d_global_t *glob,
                                     fclaw2d_patch_t *fine_patches,
                                     int blockno,
                                     int patchno)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);

    int mx,my,mz,mbc,meqn;
    double xlower,ylower,zlower,dx,dy,dz;
    double *q[4];
    int tag_patch,igrid;
    double coarsen_threshold;
    fclaw2d_patch_t *patch0;

    patch0 = &fine_patches[0];

    coarsen_threshold = gparms->coarsen_threshold;

    fclaw_clawpatch3_grid_data(glob,patch0,&mx,&my,&mz,&mbc,
                               &xlower,&ylower,&zlower,&dx,&dy,&dz);

    for (igrid = 0; igrid < 4; igrid++)
    {
        fclaw_clawpatch3_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
    }

    tag_patch = 0;
    /* 2.5 check needed */
    fclaw_clawpatch3_vt()->fort_tag4coarsening(&mx,&my,&mz,&mbc,&meqn,&xlower,&ylower,&zlower,
                                                 &dx,&dy,&dz,
                                                 &blockno, q[0],q[1],q[2],q[3],
                                                 &coarsen_threshold,&tag_patch);
    return tag_patch == 1;
}


/* -----------------------------------------------------------------
   Callback routine for tagging
   ----------------------------------------------------------------- */

void fclaw_clawpatch3_interpolate2fine(fclaw2d_global_t* glob,
                                        fclaw2d_patch_t *coarse_patch,
                                        fclaw2d_patch_t* fine_patches,
                                        int this_blockno, int coarse_patchno,
                                        int fine0_patchno)
{
    int mx,my,mz,mbc,meqn;
    double *qcoarse,*qfine;
    double *areacoarse,*areafine;
    // double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid;

    const amr_options_t* gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    fclaw2d_patch_t* fine_patch;

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

#if 0
    fclaw_clawpatch3_metric_data(glob,coarse_patch,&xp,&yp,&zp,
                                   &xd,&yd,&zd,&areacoarse);
#endif
    areacoarse = NULL;
    areafine = NULL;
    fclaw_clawpatch3_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    /* Loop over four siblings (z-ordering) */
    /* 2.5 check needed */
    for (igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw_clawpatch3_soln_data(glob,fine_patch,&qfine,&meqn);
#if 0
        if (gparms->manifold)
        {
            fclaw_clawpatch3_metric_data(glob,fine_patch,&xp,&yp,&zp,
                                          &xd,&yd,&zd,&areafine);
        }
#endif
        fclaw_clawpatch3_vt()->fort_interpolate2fine(&mx,&my,&mz,&mbc,&meqn,qcoarse,qfine,
                                                       areacoarse, areafine, &igrid,
                                                       &gparms->manifold);

    }
}

/* It seems that there is only one way to average the solution, but
   what about other things the user might want to do?   Maybe we need
   something like "average from fine" routine which handles more generic
   things, including area averaging, and maybe something to do with averaging
   stuff in aux arrays. */
void fclaw_clawpatch3_average2coarse(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *fine_patches,
                                      fclaw2d_patch_t *coarse_patch,
                                      int blockno, int fine0_patchno,
                                      int coarse_patchno)

{
    const amr_options_t* gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);
    
    int mx,my,mz,mbc, meqn;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;    
    //double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid;
    fclaw2d_patch_t *fine_patch;

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

    areacoarse = NULL;
    areafine = NULL;
#if 0
    fclaw_clawpatch3_metric_data(glob,coarse_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&areacoarse);
#endif
    fclaw_clawpatch3_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    /* 2.5 check needed */
    for(igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw_clawpatch3_soln_data(glob,fine_patch,&qfine,&meqn);
#if 0
        if (gparms->manifold)
        {
            fclaw_clawpatch3_metric_data(glob,fine_patch,&xp,&yp,&zp,
                                          &xd,&yd,&zd,&areafine);
        }
#endif
        fclaw_clawpatch3_vt()->fort_average2coarse(&mx,&my,&mz,&mbc,&meqn,qcoarse,qfine,
                                                     areacoarse, areafine, &igrid,
                                                     &gparms->manifold);

    }
}
