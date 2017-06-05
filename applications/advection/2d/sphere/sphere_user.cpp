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

#include "sphere_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <clawpack46_user_fort.h>  /* Headers for user defined fortran files */

#include <fc2d_clawpack5.h>
#include <clawpack5_user_fort.h>

#include "../all/clawpack_user.h"

void sphere_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    patch_vt->setup      = &sphere_patch_setup;  

    fclaw2d_vtable_t *vt = fclaw2d_vt();
    vt->problem_setup    = &sphere_problem_setup;  /* Version-independent */


    const user_options_t* user = sphere_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();

        clawpack46_vt->qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->rpn2      = &CLAWPACK46_RPN2ADV_MANIFOLD;
        clawpack46_vt->rpt2      = &CLAWPACK46_RPT2ADV_MANIFOLD;

        clawpatch_vt->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *clawpack5_vt = fc2d_clawpack5_vt();

        clawpack5_vt->qinit     = &CLAWPACK5_QINIT;
        clawpack5_vt->rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        clawpack5_vt->rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;

        clawpatch_vt->fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;

    }
}



#if 0
void sphere_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t *user =  sphere_user_get_options(domain);

    fclaw2d_init_vtable(&fclaw2d_vt);
    fclaw2d_vt.problem_setup = &sphere_problem_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt, &classic_claw46);

        fclaw2d_vt.patch_setup   = &sphere_patch_setup;

        /* Needed to avoid triggering refinement around block corners */
        fclaw2d_vt.fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
        fclaw2d_vt.fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

        classic_claw46.qinit = &CLAWPACK46_QINIT;
        classic_claw46.rpn2  = &CLAWPACK46_RPN2ADV_MANIFOLD;
        classic_claw46.rpt2  = &CLAWPACK46_RPT2ADV_MANIFOLD;

        fc2d_clawpack46_set_vtable(classic_claw46);

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt, &classic_claw5);

        fclaw2d_vt.patch_setup    = &sphere_patch_setup;

        /* Avoid triggering refinement around block corners */
        fclaw2d_vt.fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;
        fclaw2d_vt.fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;

        classic_claw5.qinit = &CLAWPACK5_QINIT;
        classic_claw5.rpn2  = &CLAWPACK5_RPN2ADV_MANIFOLD;
        classic_claw5.rpt2  = &CLAWPACK5_RPT2ADV_MANIFOLD;

        fc2d_clawpack5_set_vtable(classic_claw5);
    }
    fclaw2d_set_vtable(domain,&fclaw2d_vt);
}
#endif


void sphere_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t *user = sphere_get_options(glob);

    SPHERE_SETPROB(&user->revs_per_second);
}


void sphere_patch_setup(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    const user_options_t* user = sphere_get_options(glob);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    if (user->claw_version == 4)
    {
        USER46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                               &maux,aux,&this_block_idx,xd,yd,zd,area);
    }
    else if(user->claw_version == 5)
    {
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                              &maux,aux,&this_block_idx,xd,yd,zd,area);
    }
}
