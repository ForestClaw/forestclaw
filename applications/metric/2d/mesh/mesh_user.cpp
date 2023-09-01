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

#include "mesh_user.h"
#include <fclaw_forestclaw.h>
#include <fclaw_clawpatch.h>


/* Two versions of Clawpack */
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

void mesh_link_solvers(fclaw_global_t *glob)
{
    const user_options_t* user = mesh_user_get_options(glob);

    fclaw_vt(glob)->problem_setup = &mesh_problem_setup;  /* Version-independent */

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vt(glob)->qinit     = &CLAWPACK46_QINIT;
        fc2d_clawpack46_vt(glob)->setaux    = &CLAWPACK46_SETAUX;
        fc2d_clawpack46_vt(glob)->rpn2      = &CLAWPACK46_RPN2ADV;
        fc2d_clawpack46_vt(glob)->rpt2      = &CLAWPACK46_RPT2ADV;
        fc2d_clawpack46_vt(glob)->b4step2   = &CLAWPACK46_B4STEP2;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vt(glob)->qinit     = &CLAWPACK5_QINIT;
        fc2d_clawpack5_vt(glob)->setaux    = &CLAWPACK5_SETAUX;
        fc2d_clawpack5_vt(glob)->b4step2   = &CLAWPACK5_B4STEP2;
        fc2d_clawpack5_vt(glob)->rpn2      = &CLAWPACK5_RPN2ADV;
        fc2d_clawpack5_vt(glob)->rpt2      = &CLAWPACK5_RPT2ADV;
    }
}

void mesh_problem_setup(fclaw_global_t* glob)
{
    const user_options_t* user = mesh_user_get_options(glob);

    double period = user->period;
    MESH_SETPROB(&period);
}
