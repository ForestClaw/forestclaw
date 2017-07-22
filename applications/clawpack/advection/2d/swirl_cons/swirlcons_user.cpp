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

#include "swirlcons_user.h"

#include <fclaw2d_include_all.h>

/* Two versions of Clawpack */
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>  /* Headers for user defined fortran files */

#include "../all/advection_user_fort.h"

void swirlcons_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    fc2d_clawpack46_options_t *clawopt = fc2d_clawpack46_get_options(glob);

    vt->problem_setup = &swirlcons_problem_setup;  /* Version-independent */

    const user_options_t* user = swirlcons_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();
        clawpack46_vt->qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->setaux    = &CLAWPACK46_SETAUX;
        if (user->cons_rp)
        {
            clawopt->use_fwaves = 0;
            clawpack46_vt->rpn2      = &RPN2CONS_CC;                                
        }
        else if (clawopt->use_fwaves)
        {
            clawpack46_vt->rpn2      = &RPN2FWAVE;                
        }
        else
        {
            clawopt->use_fwaves = 0;
            clawpack46_vt->rpn2      = &RPN2CONS;                                
        }
        clawpack46_vt->rpt2      = &RPT2CONS_CC;
        clawpack46_vt->b4step2   = &SWIRLCONS_B4STEP2;
    }
    else if (user->claw_version == 5)
    {
        fclaw_global_essentialf("swirl (cons) : Clawpack 5.x implementation not implemented\n");
        exit(0);
    }
}

void swirlcons_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = swirlcons_get_options(glob);

    double period = user->period;
    int ex = user->example;
    SWIRL_SETPROB(&period,&ex);
}





