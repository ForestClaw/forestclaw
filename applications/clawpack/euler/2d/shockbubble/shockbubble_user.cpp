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

#include "shockbubble_user.h"

#include <fclaw_include_all.h>

#include <fclaw_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>

#include <fc2d_clawpack5.h>
#include <fc2d_clawpack5_options.h>

#include "../rp/euler_user_fort.h"


void shockbubble_problem_setup(fclaw_global_t* glob)
{
    const user_options_t* user = shockbubble_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",   user->idisc,"\% idisc\n");
        fprintf(f,  "%-24.16f   %s",   user->gamma,"\% gamma\n");
        fprintf(f,  "%-24.16f   %s",user->x0,"\% x0\n");
        fprintf(f,  "%-24.16f   %s",user->y0,"\% y0\n");
        fprintf(f,  "%-24.16f   %s",user->r0,"\% r0\n");
        fprintf(f,  "%-24.16f   %s",user->rhoin,"\% rhoin\n");
        fprintf(f,  "%-24.16f   %s",user->pinf,"\% pinf\n");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
    fclaw_domain_barrier (glob->domain);  /* redundant?  */
    SETPROB();
}

void shockbubble_link_solvers(fclaw_global_t *glob)
{
    const user_options_t* user = shockbubble_get_options(glob);
    fclaw_vtable_t *fc_vt = fclaw_vt(glob);

    fc_vt->problem_setup = &shockbubble_problem_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);
        fc2d_clawpack46_options_t *clawopt = fc2d_clawpack46_get_options(glob);

        claw46_vt->fort_qinit  = &CLAWPACK46_QINIT;
        claw46_vt->fort_setaux = &CLAWPACK46_SETAUX;
        claw46_vt->fort_bc2    = &CLAWPACK46_BC2;   /* Special  BCs at left edge */
        claw46_vt->fort_src2   = &CLAWPACK46_SRC2;  /* To simulate axis-symmetric */

        switch (clawopt->mwaves)
        {
        case 4:
            /* Requires meqn=4 */
            claw46_vt->fort_rpn2   = &CLAWPACK46_RPN2_EULER4;  /* No tracer */
            claw46_vt->fort_rpt2   = &CLAWPACK46_RPT2_EULER4;
            break;
        case 5:
            /* Requires meqn=5 */
            claw46_vt->fort_rpn2   = &CLAWPACK46_RPN2_EULER5;  /* Includes a tracer */
            claw46_vt->fort_rpt2   = &CLAWPACK46_RPT2_EULER5;
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt(glob);
        fc2d_clawpack5_options_t *clawopt = fc2d_clawpack5_get_options(glob);

        claw5_vt->fort_qinit  = &CLAWPACK5_QINIT;
        claw5_vt->fort_setaux = &CLAWPACK5_SETAUX;
        claw5_vt->fort_bc2    = &CLAWPACK5_BC2;   /* Special  BCs at left edge */
        claw5_vt->fort_src2   = &CLAWPACK5_SRC2;  /* To simulate axis-symmetric */
      switch (clawopt->mwaves)
        {
        case 4:
            /* Requires meqn=4 */
            claw5_vt->fort_rpn2   = &CLAWPACK5_RPN2_EULER4;  /* no tracer */
            claw5_vt->fort_rpt2   = &CLAWPACK5_RPT2_EULER4;
            break;
        case 5:
            /* Requires meqn=5 */
            claw5_vt->fort_rpn2   = &CLAWPACK5_RPN2_EULER5;  /* Includes a tracer */
            claw5_vt->fort_rpt2   = &CLAWPACK5_RPT2_EULER5;
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }
    }
}