/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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


#include "replicated_user.h"

#if 0
#include <fclaw2d_include_all.h>

/* Two versions of Clawpack */
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>
#endif

#include "../all/advection_user.h"

void replicated_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();

    vt->problem_setup = &replicated_problem_setup; 

    const user_options_t* user = replicated_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();
        claw46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        claw46_vt->fort_setaux    = &CLAWPACK46_SETAUX;
        claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
        claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;

#if 0
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        clawpatch_vt->fort_tag4coarsening = &CLAWPATCH46_TAG4COARSENING;
        clawpatch_vt->fort_tag4refinement = &CLAWPATCH46_TAG4REFINEMENT;
#endif        
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt();
        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        claw5_vt->fort_setaux    = &CLAWPACK5_SETAUX;
        claw5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV;
        claw5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV;

#if 0
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        clawpatch_vt->fort_tag4coarsening = &CLAWPATCH5_TAG4COARSENING;
        clawpatch_vt->fort_tag4refinement = &CLAWPATCH5_TAG4REFINEMENT;        
#endif        
    }
}

void replicated_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = replicated_get_options(glob);

    int example = user->example;  /* Macros don't expand properly without this */
    REPLICATED_SETPROB(&example);    
}


