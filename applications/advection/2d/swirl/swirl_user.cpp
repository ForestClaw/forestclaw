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

#include "swirl_user.h"
#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>
#include <fc2d_dummy.h>

static fclaw2d_vtable_t fclaw2d_vt;
static fc2d_clawpack46_vtable_t classic_claw46;
static fc2d_clawpack5_vtable_t classic_claw5;

void swirl_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw_app_t* app;
    user_options_t* user;

    app = fclaw2d_domain_get_app(domain);
    user = (user_options_t*) fclaw_app_get_user(app);

    fclaw2d_init_vtable(&fclaw2d_vt);

    fclaw2d_vt.problem_setup            = &swirl_problem_setup;
    fclaw2d_vt.patch_setup              = &swirl_patch_setup;    /* Needs to call setaux */

    if (user->claw_version == 4)
    {
        /* Needed for the clawpack46 package */
        fc2d_clawpack46_init_vtable(&fclaw2d_vt,&classic_claw46);

        /* Customized refinement so that initial conditions are properly tagged. */
        fclaw2d_vt.fort_tag4refinement      = &CLAWPACK46_TAG4REFINEMENT;

        classic_claw46.qinit     = &CLAWPACK46_QINIT;
        classic_claw46.setaux    = &CLAWPACK46_SETAUX;
        classic_claw46.rpn2      = &CLAWPACK46_RPN2;
        classic_claw46.rpt2      = &CLAWPACK46_RPT2;
        classic_claw46.b4step2   = &CLAWPACK46_B4STEP2;

        fc2d_clawpack46_set_vtable(&classic_claw46);

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_init_vtable(&fclaw2d_vt,&classic_claw5);

        /* Customized refinement so that initial conditions are properly tagged. */
        fclaw2d_vt.fort_tag4refinement   = &CLAWPACK5_TAG4REFINEMENT;

        classic_claw5.qinit     = &CLAWPACK5_QINIT;
        classic_claw5.setaux    = &CLAWPACK5_SETAUX;
        classic_claw5.b4step2   = &CLAWPACK5_B4STEP2;
        classic_claw5.rpn2      = &CLAWPACK5_RPN2;
        classic_claw5.rpt2      = &CLAWPACK5_RPT2;

        fc2d_clawpack5_set_vtable(&classic_claw5);
    }

    fclaw2d_set_vtable(domain,&fclaw2d_vt);
}


void swirl_problem_setup(fclaw2d_domain_t* domain)
{
    fclaw_app_t* app;
    user_options_t* user;

    app = fclaw2d_domain_get_app(domain);
    user = (user_options_t*) fclaw_app_get_user(app);

    SWIRL_SETPROB(&user->period);
}

void swirl_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    fclaw_app_t* app;
    user_options_t* user;

    app = fclaw2d_domain_get_app(domain);
    user = (user_options_t*) fclaw_app_get_user(app);

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        return;
    }
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_setaux(domain,this_patch,this_block_idx,this_patch_idx);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_setaux(domain,this_patch,this_block_idx,this_patch_idx);
    }
    /* Dummy setup - use multiple libraries */
    fc2d_dummy_setup_patch(domain,this_patch,this_block_idx,this_patch_idx);
}

int swirl_patch_tag4refinement(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int blockno, int this_patch_idx,
                               int initflag)
{
    fclaw2d_vtable_t vt;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;
    const amr_options_t *amropt;
    double rt;

    amropt = get_domain_parms(domain);
    rt = amropt->refine_threshold;

    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    tag_patch = 0;
    vt.fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,
                           &dx,&dy,&blockno, q,&rt,&initflag,
                           &tag_patch);
    return tag_patch;
}
