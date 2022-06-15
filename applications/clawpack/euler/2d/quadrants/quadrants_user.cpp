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

#include "quadrants_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

#include "../rp/euler_user_fort.h"

static
void quadrants_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = quadrants_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24.16f   %s",   user->gamma,"\% gamma\n");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
    fclaw2d_domain_barrier (glob->domain);  /* redundant?  */
    SETPROB();
}

void quadrants_link_solvers(fclaw2d_global_t *glob)
{
    user_options_t *user =  quadrants_get_options(glob);
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt(glob);

    fclaw_vt->problem_setup = &quadrants_problem_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);

        claw46_vt->fort_qinit = &CLAWPACK46_QINIT;
        claw46_vt->fort_rpn2  = &CLAWPACK46_RPN2_EULER4;
        claw46_vt->fort_rpt2  = &CLAWPACK46_RPT2_EULER4;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt(glob);

        claw5_vt->fort_qinit = &CLAWPACK5_QINIT;
        claw5_vt->fort_rpn2  = &CLAWPACK5_RPN2_EULER4;  /* Signature is unchanged */
        claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2_EULER4;
    }
}

