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

#include "radial_user.h"
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

static fclaw2d_vtable_t fclaw2d_vt;
static fc2d_clawpack46_vtable_t classic_claw46;
static fc2d_clawpack5_vtable_t classic_claw5;


void radial_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t* user = radial_user_get_options(domain);

    fclaw2d_init_vtable(&fclaw2d_vt);
    fclaw2d_vt.problem_setup = &radial_problem_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt, &classic_claw46);
        fclaw2d_vt.patch_setup = &radial_patch_setup;

        classic_claw46.qinit = &CLAWPACK46_QINIT;
        if (user->example == 0)
        {
            classic_claw46.rpn2  = &CLAWPACK46_RPN2;
            classic_claw46.rpt2  = &CLAWPACK46_RPT2;
        }
        else if (user->example == 1)
        {
            fclaw_global_essentialf("Not implemented yet\n");
            exit(0);
        }

        fclaw2d_vt.fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;  /* User defined */
        fclaw2d_vt.fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;  /* User defined */

        fc2d_clawpack46_set_vtable(classic_claw46);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt, &classic_claw5);
        fclaw2d_vt.patch_setup = &radial_patch_setup;

        classic_claw5.qinit = &CLAWPACK5_QINIT;
        if (user->example == 0)
        {
            classic_claw5.rpn2  = &CLAWPACK5_RPN2;
            classic_claw5.rpt2  = &CLAWPACK5_RPT2;
        }
        else if (user->example == 1)
        {
            classic_claw5.rpn2  = &CLAWPACK5_RPN2_MANIFOLD;
            classic_claw5.rpt2  = &CLAWPACK5_RPT2_MANIFOLD;
        }

        fclaw2d_vt.fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;  /* User defined */
        fclaw2d_vt.fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;  /* User defined */

        fc2d_clawpack5_set_vtable(classic_claw5);
    }

    fclaw2d_set_vtable(domain,&fclaw2d_vt);
}

void radial_problem_setup(fclaw2d_domain_t* domain)
{
    user_options_t* user = radial_user_get_options(domain);

    /* rho, bulk are inputs; cc and zz are outputs.  Results are
       stored in a common block */
    RADIAL_SETPROB(&user->rho,&user->bulk,&user->cc,&user->zz);
}

void radial_patch_setup(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    const user_options_t* user = radial_user_get_options(domain);

    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *xnormals,*ynormals,*xtangents,*ytangents;
    double *surfnormals,*edgelengths,*curvature;

    if (user->example == 0)
    {
        /* Not mapped */
        return;
    }

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        /* Mapped info is needed only for an update */
        return;
    }

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_data2(domain,this_patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,&edgelengths,
                                   &curvature);

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_define_auxarray(domain,this_patch);
        fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);
        USER46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                               &dx,&dy,&maux,aux,
                               xnormals,ynormals,edgelengths,area);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_define_auxarray(domain,this_patch);
        fc2d_clawpack5_aux_data(domain,this_patch,&aux,&maux);
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                              &dx,&dy,&maux,aux,
                              xnormals,ynormals,edgelengths,area);
    }
}
