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
    fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();

    clawpack46_vt->qinit     = &CLAWPACK46_QINIT;
    clawpack46_vt->setaux    = &CLAWPACK46_SETAUX;
    clawpack46_vt->rpt2      = &RPT2CONS;

    switch(user->rp_solver)
    {
        case 1:
            clawopt->use_fwaves = 0;
            clawpack46_vt->rpn2      = &RPN2CONS_QS; 
            break; 

        case 2:
            clawopt->use_fwaves = 0;
            clawpack46_vt->rpn2      = &RPN2CONS_WD; 
            break; 

        case 3:
            clawopt->use_fwaves = 0;
            clawpack46_vt->rpn2      = &RPN2CONS_EC; 
            break;

        case 4:
            clawopt->use_fwaves = 1;
            clawpack46_vt->rpn2      = &RPN2CONS_FW; 
            break;
    }
 }

void swirlcons_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = swirlcons_get_options(glob);

    int ex = user->example;
    SWIRL_SETPROB(&ex);
}





