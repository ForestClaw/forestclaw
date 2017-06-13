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

#include "interface_user.h"

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <clawpack46_user_fort.h>  /* Headers for user defined fortran files */

#include <fc2d_clawpack5.h>
#include <clawpack5_user_fort.h>

#include "../rp/clawpack_user.h"

void interface_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    fclaw2d_vtable_t *vt = fclaw2d_vt();

    vt->problem_setup = &interface_problem_setup;  /* Version-independent */

    const user_options_t* user = interface_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();

        clawpack46_vt->qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->setaux    = &CLAWPACK46_SETAUX;
        clawpack46_vt->rpn2      = &CLAWPACK46_RPN2;
        clawpack46_vt->rpt2      = &CLAWPACK46_RPT2;

        clawpatch_vt->fort_tag4refinement  = &CLAWPACK46_TAG4REFINEMENT;  

        /* Use default coarsening: coarsening only if solution is sufficiently flat */
        clawpatch_vt->fort_tag4coarsening  = &CLAWPACK46_TAG4COARSENING;  

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t    *clawpack5_vt = fc2d_clawpack5_vt();

        clawpack5_vt->qinit   = &CLAWPACK5_QINIT;
        clawpack5_vt->setaux  = &CLAWPACK5_SETAUX;
        clawpack5_vt->rpn2    = &CLAWPACK5_RPN2;
        clawpack5_vt->rpt2    = &CLAWPACK5_RPT2;

        clawpatch_vt->fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;

        /* Use default coarsening: coarsening only if solution is sufficiently flat */
        clawpatch_vt->fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;
    }
}


#if 0
void interface_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t* user = interface_user_get_options(domain);

    fclaw2d_init_vtable(&fclaw2d_vt);
    fclaw2d_vt.problem_setup = &interface_setup_problem;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt, &classic_claw46);

        classic_claw46.qinit      = &CLAWPACK46_QINIT;
        classic_claw46.setaux     = &CLAWPACK46_SETAUX;
        classic_claw46.rpn2       = &CLAWPACK46_RPN2;
        classic_claw46.rpt2       = &CLAWPACK46_RPT2;

        fclaw2d_vt.fort_tag4refinement  = &CLAWPACK46_TAG4REFINEMENT;  /* User defined */

        fc2d_clawpack46_set_vtable(classic_claw46);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt, &classic_claw5);

        classic_claw5.qinit      = &CLAWPACK5_QINIT;
        classic_claw5.setaux     = &CLAWPACK5_SETAUX;
        classic_claw5.rpn2       = &CLAWPACK5_RPN2;
        classic_claw5.rpt2       = &CLAWPACK5_RPT2;

        fclaw2d_vt.fort_tag4refinement      = &CLAWPACK5_TAG4REFINEMENT;  /* User defined */

        fc2d_clawpack5_set_vtable(classic_claw5);
    }

    fclaw2d_set_vtable(domain,&fclaw2d_vt);
}
#endif

void interface_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = interface_get_options(glob);
    INTERFACE_SETPROB(&user->rhol,&user->cl,&user->rhor,&user->cr);
}
