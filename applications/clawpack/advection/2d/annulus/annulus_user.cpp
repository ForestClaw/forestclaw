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

#include "annulus_user.h"

static
void annulus_problem_setup(fclaw_global_t *glob)
{
    const user_options_t *user = annulus_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24.16f   %s",user->beta,"\% beta\n");
        fclose(f);
    }

    /* Make sure node 0 writes 'setprob.data' before proceeding */
    fclaw_domain_barrier (glob->domain);  /* redundant?  */
 
    SETPROB();
}


static
void annulus_patch_setup(fclaw_global_t *glob,
                         fclaw_patch_t *patch,
                         int blockno,
                         int patchno)
{
    const user_options_t* user = annulus_get_options(glob);
    advection_patch_setup_manifold(glob,patch,blockno,patchno,
                                   user->claw_version);
}


void annulus_link_solvers(fclaw_global_t *glob)
{
    const user_options_t   *user      =  annulus_get_options(glob);

    fclaw_vtable_t           *vt           = fclaw_vt(glob);
    fclaw_patch_vtable_t     *patch_vt     = fclaw_patch_vt(glob);

    vt->problem_setup  = &annulus_problem_setup;
    patch_vt->setup    = &annulus_patch_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt(glob);
        clawpack46_vt->fort_qinit   = CLAWPACK46_QINIT;
        clawpack46_vt->fort_rpn2    = CLAWPACK46_RPN2ADV_MANIFOLD;
        clawpack46_vt->fort_rpt2    = CLAWPACK46_RPT2ADV_MANIFOLD;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt(glob);
        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        claw5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        claw5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;
    }
}


