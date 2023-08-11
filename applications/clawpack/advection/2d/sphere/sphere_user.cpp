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

#include "sphere_user.h"

static
void sphere_problem_setup(fclaw_global_t* glob)
{
    const user_options_t *user = sphere_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->revs_per_second,"\% revs_per_second");
        fclose(f);
    }

    fclaw2d_domain_barrier(glob->domain);
    SETPROB();
}

static
void sphere_patch_setup(fclaw_global_t *glob,
                        fclaw_patch_t *patch,
                        int blockno,
                        int patchno)
{
    const user_options_t* user = sphere_get_options(glob);
    advection_patch_setup_manifold(glob,patch,blockno,patchno,
                                   user->claw_version);
}


void sphere_link_solvers(fclaw_global_t *glob)
{
    fclaw_vtable_t *vt = fclaw_vt(glob);
    vt->problem_setup    = &sphere_problem_setup;  /* Version-independent */

    fclaw_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
    patch_vt->setup      = &sphere_patch_setup;  

    const user_options_t* user = sphere_get_options(glob);
    if (user->example == 1)
    {
        /* Needed to get correct handling of block corners */
        fclaw_clawpatch_use_pillowsphere(glob);
    }


    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);

        claw46_vt->fort_qinit  = CLAWPACK46_QINIT;
        claw46_vt->fort_rpn2   = CLAWPACK46_RPN2ADV_MANIFOLD;
        claw46_vt->fort_rpt2   = CLAWPACK46_RPT2ADV_MANIFOLD;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt(glob);

        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        claw5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        claw5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;
    }
}
