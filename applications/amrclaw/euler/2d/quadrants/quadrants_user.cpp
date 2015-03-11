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

#include "quadrants_user.H"
#include "fclaw2d_clawpatch.H"
#include "fc2d_clawpack46.H"


static fclaw2d_vtable_t vt;
static fc2d_clawpack46_vtable_t classic_claw;

void quadrants_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);
    fc2d_clawpack46_init_vtable(&classic_claw);

    vt.problem_setup = &quadrants_problem_setup;
    /* Don't explicitly set a "setprob" function unless it has the same
       signature as the default (i.e. no arguments).  */
    /* classic_claw.setprob = &SETPROB; */

    vt.patch_initialize = &fc2d_clawpack46_qinit;
    classic_claw.qinit = &QINIT;

    vt.patch_physical_bc = &fc2d_clawpack46_bc2;  /* Set to bc2 by default */

    vt.patch_single_step_update = &fc2d_clawpack46_update;
    classic_claw.rpn2 = &RPN2EU3;  /* Signature is unchanged */
    classic_claw.rpt2 = &RPT2;

    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);
}

void quadrants_problem_setup(fclaw2d_domain_t* domain)
{
    const user_options_t* user;
    user = (user_options_t*) fclaw2d_domain_get_user_options(domain);

    QUADRANTS_SETPROB(&user->gamma);
}
