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

#include "bump_user.h"

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

#include "../rp/shallow_user_fort.h"

void bump_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();

    vt->problem_setup = &bump_problem_setup;  /* Version-independent */

    const user_options_t* user = bump_get_options(glob);
    if (user->claw_version == 4)
    {
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

        claw46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2;
        claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2;

        /* Avoid tagging block corners in 5 patch example*/
        clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;


    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t    *claw5_vt = fc2d_clawpack5_vt();

        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;

        if (user->example == 0)
        {
            claw5_vt->fort_rpn2 = &CLAWPACK5_RPN2;
            claw5_vt->fort_rpt2 = &CLAWPACK5_RPT2;
        }
        else if (user->example == 1)
        {
            fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
            fclaw2d_patch_vtable_t         *patch_vt = fclaw2d_patch_vt();

            patch_vt->setup = &bump_patch_setup;

            claw5_vt->fort_rpn2  = &CLAWPACK5_RPN2_MANIFOLD;
            claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2_MANIFOLD;

            /* Avoid tagging block corners in 5 patch example*/
            clawpatch_vt->fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;
            clawpatch_vt->fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;
        }
    }
}


void bump_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = bump_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24.16f   %s",user->gravity,"\% gravity\n");
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    BUMP_SETPROB();
}


void bump_patch_setup(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *xnormals,*ynormals,*xtangents,*ytangents;
    double *surfnormals,*edgelengths,*curvature;

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        /* Mapped info is needed only for an update */
        return;
    }

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_data2(glob,this_patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,&edgelengths,
                                   &curvature);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    
    USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                          &dx,&dy,&maux,aux,
                          xnormals,xtangents,
                          ynormals,ytangents,
                          surfnormals,area);
}
