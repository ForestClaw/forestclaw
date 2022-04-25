/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

static
void swirl_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = swirl_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->period,"\% period");
        fclose(f);
    }

    /* Make sure node 0 has written 'setprob.data' before proceeding */
    fclaw2d_domain_barrier (glob->domain);

    SETPROB();
}

void swirl_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    vt->problem_setup = &swirl_problem_setup;  /* Version-independent */

    swirl_initialize_rays(glob);

    const user_options_t* user = swirl_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();        

        clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
        clawpack46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;

        /* Velocity is set here rather than in setaux, because we have a 
           time dependent velocity field */
        clawpack46_vt->fort_b4step2   = &CLAWPACK46_B4STEP2;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *clawpack5_vt = fc2d_clawpack5_vt();

        clawpack5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        clawpack5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV;
        clawpack5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV;  

        /* Velocity is set here rather than in setaux, because we have a 
           time dependent velocity field */
        clawpack5_vt->fort_b4step2   = &CLAWPACK5_B4STEP2;
    }
}





