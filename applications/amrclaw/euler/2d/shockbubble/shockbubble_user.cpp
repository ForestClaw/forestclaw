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

#include "shockbubble_user.h"

#include <fclaw2d_clawpatch.h>
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

static fclaw2d_vtable_t fclaw2d_vt;
static fc2d_clawpack46_vtable_t classic_claw46;
static fc2d_clawpack5_vtable_t classic_claw5;

void shockbubble_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t* user = shockbubble_user_get_options(domain);

    fclaw2d_init_vtable(&fclaw2d_vt);
    fclaw2d_vt.problem_setup = &shockbubble_problem_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt, &classic_claw46);
        fc2d_clawpack46_options_t *clawopt = fc2d_clawpack46_get_options(domain);

#if 0
        fclaw2d_vt.patch_setup = &shockbubble_patch_setup;
#endif

        classic_claw46.qinit  = &CLAWPACK46_QINIT;
        classic_claw46.setaux = &CLAWPACK46_SETAUX;
        classic_claw46.bc2    = &CLAWPACK46_BC2;  /* Special input BCs */
        classic_claw46.src2   = &CLAWPACK46_SRC2;  /* To simulate axis-symmetric */

        switch (clawopt->mwaves)
        {
        case 3:
        case 4:
            fclaw_global_essentialf("3-wave and 4-wave solvers not implemented for Clawpack 4.6\n");
            exit(0);
            break;
        case 5:
            classic_claw46.rpn2   = &CLAWPACK46_RPN2_EULER5;  /* Includes a tracer */
            classic_claw46.rpt2   = &CLAWPACK46_RPT2_EULER5;  /* Includes a tracer */
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* Use divided differences to tag grids */
        fclaw2d_vt.fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
        fclaw2d_vt.fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

        fc2d_clawpack46_set_vtable(classic_claw46);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt, &classic_claw5);
        fc2d_clawpack5_options_t *clawopt = fc2d_clawpack5_get_options(domain);

        classic_claw5.qinit  = &CLAWPACK5_QINIT;
        classic_claw5.setaux = &CLAWPACK5_SETAUX;
        classic_claw5.bc2    = &CLAWPACK5_BC2;  /* Special input BCs */
        classic_claw5.src2   = &CLAWPACK5_SRC2;  /* To simulate axis-symmetric */

        switch (clawopt->mwaves)
        {
        case 3:
        case 4:
            fclaw_global_essentialf("3-wave and 4-wave solvers not implemented for Clawpack 5.0\n");
            exit(0);
            break;
        case 5:
            classic_claw5.rpn2   = &CLAWPACK5_RPN2_EULER5;  /* Includes a tracer */
            classic_claw5.rpt2   = &CLAWPACK5_RPT2_EULER5;
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* Use divided differences to tag grids */
        fclaw2d_vt.fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;
        fclaw2d_vt.fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;

        fc2d_clawpack5_set_vtable(classic_claw5);
    }
    fclaw2d_set_vtable(domain,&fclaw2d_vt);

#if 0
    classic_claw.src2 = &SRC2;
    classic_claw.rpn2 = &RPN2EU5;  /* Signature is unchanged */
    classic_claw.rpt2 = &RPT2EU5;
#endif
}

void shockbubble_problem_setup(fclaw2d_domain_t* domain)
{
    const user_options_t* user = shockbubble_user_get_options(domain);

    SHOCKBUBBLE_SETPROB(&user->gamma, &user->x0, &user->y0, &user->r0,
                        &user->rhoin, &user->pinf);
}

#if 0
void shockbubble_patch_setup(fclaw2d_domain_t* domain,
                             fclaw2d_patch_t* this_patch,
                             int blockno,
                             int patchno)
{
    /* I don't have the scaling right on this problem, and just setting
       mcapa doesn't really fix the issue.   I need a way to set [ax,ay,bx,by]
       inside of ClawPatch */
    fc2d_clawpack46_setaux(domain,this_patch,blockno,patchno);
    fc2d_clawpack46_set_capacity(domain,this_patch,blockno,patchno);
}
#endif
