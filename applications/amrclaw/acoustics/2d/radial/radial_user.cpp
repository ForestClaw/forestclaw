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
    fclaw2d_vt.problem_setup = &radial_setup_problem;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt, &classic_claw46);

        classic_claw46.qinit = &CLAWPACK46_QINIT;
        classic_claw46.rpn2  = &CLAWPACK46_RPN2;
        classic_claw46.rpt2  = &CLAWPACK46_RPT2;


        fclaw2d_vt.fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;  /* User defined */

        fc2d_clawpack46_set_vtable(classic_claw46);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt, &classic_claw5);

        classic_claw5.qinit = &CLAWPACK5_QINIT;
        classic_claw5.rpn2  = &CLAWPACK5_RPN2;
        classic_claw5.rpt2  = &CLAWPACK5_RPT2;

        fclaw2d_vt.fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;  /* User defined */

        fc2d_clawpack5_set_vtable(classic_claw5);
    }

    fclaw2d_set_vtable(domain,&fclaw2d_vt);

#if 0
    fclaw2d_init_vtable(&fclaw2d_vt);


    fc2d_clawpack46_init_vtable(&classic_claw);

    vt.problem_setup = &radial_problem_setup;  /* setprob called from here */
    /* classic_claw.setprob = &SETPROB; */

    vt.patch_initialize = &fc2d_clawpack46_qinit;
    classic_claw.qinit = &QINIT;

    vt.patch_physical_bc = fc2d_clawpack46_bc2;

    vt.patch_single_step_update = &fc2d_clawpack46_update;
    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;

    /* Refine if a patch contain values that exceed (in magnitude)
       a refinement threshold (in fclaw_options.ini) */
    vt.fort_tag4refinement = &TAG4REFINEMENT;


    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);
#endif
}

void radial_setup_problem(fclaw2d_domain_t* domain)
{
    user_options_t* user = radial_user_get_options(domain);

    /* rho, bulk are inputs; cc and zz are outputs.  Results are
       stored in a common block */
    RADIAL_SETPROB(&user->rho,&user->bulk,&user->cc,&user->zz);
}
