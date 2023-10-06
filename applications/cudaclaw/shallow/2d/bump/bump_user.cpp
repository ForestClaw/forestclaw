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

#include "bump_user.h"

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>


#include "../../../../clawpack/shallow/2d/rp/shallow_user_fort.h"


static
void bump_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = bump_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        if(user->cuda != 0)
        {
            fprintf(f,  "%-24.16f   %s",user->gravity,"\% gravity\n");
        }
        else
        {
            fprintf(f,  "%-24d   %s",   user->example,"\% example\n");
            fprintf(f,  "%-24.16f   %s",user->gravity,"\% gravity\n");
        }
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);

    if(user->cuda != 0)
    {
        setprob_cuda();
    }
    else
    {
       BUMP_SETPROB(); 
    }
}


void bump_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &bump_problem_setup;  /* Version-independent */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);

    const user_options_t* user = bump_get_options(glob);
    if(user->cuda == 0)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);
        claw46_vt->fort_qinit     = &CUDACLAW_QINIT;
        claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2;
        claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2;
        claw46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE;
    }
    else
    {
        fc2d_cudaclaw_vtable_t *cuclaw_vt = fc2d_cudaclaw_vt(glob);
        cuclaw_vt->fort_qinit  = &CUDACLAW_QINIT;

        bump_assign_rpn2(&cuclaw_vt->cuda_rpn2);
        FCLAW_ASSERT(cuclaw_vt->cuda_rpn2 != NULL);

        bump_assign_rpt2(&cuclaw_vt->cuda_rpt2);
        FCLAW_ASSERT(cuclaw_vt->cuda_rpt2 != NULL);

        bump_assign_speeds(&cuclaw_vt->cuda_speeds);
        FCLAW_ASSERT(cuclaw_vt->cuda_speeds != NULL);
    }
}


#if 0
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
#endif
