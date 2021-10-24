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

#include "swirl_user.h"

void swirl_link_solvers(fclaw2d_global_t *glob)
{
    /* example of how to set up a user defined criteria */
    fclaw3dx_clawpatch_vtable_t *clawpatch_vt = fclaw3dx_clawpatch_vt();
    clawpatch_vt->fort_user_exceeds_threshold = &USER_EXCEEDS_TH;

    const user_options_t* user = swirl_get_options(glob);
    if (user->claw_version == 4)
    {
        fc3d_clawpack46_vtable_t *clawpack46_vt = fc3d_clawpack46_vt();        

        clawpack46_vt->fort_setprob   = &SETPROB;
        clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->fort_setaux    = &CLAWPACK46_SETAUX;
        clawpack46_vt->fort_rpn3      = &CLAWPACK46_RPN3;
        clawpack46_vt->fort_rpt3      = &CLAWPACK46_RPT3;
        clawpack46_vt->fort_rptt3      = &CLAWPACK46_RPTT3;
        // clawpack46_vt->fort_b4step3   = &CLAWPACK46_B4STEP3;        
    }
    else if (user->claw_version == 5)
    {
        printf("swirl_user.cpp : Example not implemented for Claw version 5.\n");
        exit(0);
    }
}





