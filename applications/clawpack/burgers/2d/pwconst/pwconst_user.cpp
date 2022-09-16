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

#include "pwconst_user.h"

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <clawpack46_user_fort.h>  /* Headers for user defined fortran files */

#include <fc2d_clawpack5.h>
#include <clawpack5_user_fort.h>

void pwconst_link_solvers(fclaw2d_global_t *glob)
{
    const user_options_t* user = pwconst_get_options(glob);
    
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);

        claw46_vt->fort_qinit  = &CLAWPACK46_QINIT;
        claw46_vt->fort_rpn2   = &CLAWPACK46_RPN2;
        claw46_vt->fort_rpt2   = &CLAWPACK46_RPT2;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t    *claw5_vt = fc2d_clawpack5_vt(glob);

        claw5_vt->fort_qinit = &CLAWPACK5_QINIT;
        claw5_vt->fort_rpn2  = &CLAWPACK5_RPN2;
        claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2;
    }
}



#if 0
void pwconst_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t* user = pwconst_user_get_options(domain);

    fclaw2d_init_vtable(&fclaw2d_vt);

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt,&classic_claw46);

        classic_claw46.qinit     = &CLAWPACK46_QINIT;
        classic_claw46.rpn2      = &CLAWPACK46_RPN2;
        classic_claw46.rpt2      = &CLAWPACK46_RPT2;

        fc2d_clawpack46_set_vtable(classic_claw46);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt,&classic_claw5);

        classic_claw5.qinit     = &CLAWPACK5_QINIT;
        classic_claw5.rpn2      = &CLAWPACK5_RPN2;
        classic_claw5.rpt2      = &CLAWPACK5_RPT2;

        fc2d_clawpack5_set_vtable(classic_claw5);
    }

    fclaw2d_set_vtable(domain,&fclaw2d_vt);

#if 0
    fclaw2d_init_vtable(&fclaw2d_vt);
    fc2d_clawpack46_init_vtable(&classic_claw);

    vt.patch_initialize           = &fc2d_clawpack46_qinit;
    classic_claw.qinit = &QINIT;

    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;

    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);
#endif
}
#endif
