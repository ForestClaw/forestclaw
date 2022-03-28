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

#include "torus_user.h"

#include "../all/advection_user.h"

static
void torus_problem_setup(fclaw2d_global_t *glob)
{
    const user_options_t* user = torus_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,"\% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->beta,"\% beta\n");
        fprintf(f,  "%-24.6f   %s",user->revs_per_s,"\% revs_per_second\n");
        fclose(f);
    }
    /* Make sure node 0 writes `setprob.data` before proceeding */
    fclaw2d_domain_barrier (glob->domain);
    
    SETPROB();
}

static
void torus_patch_setup(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *patch,
                       int blockno,
                       int patchno)
{
    const user_options_t* user = torus_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_setaux(glob,patch,blockno,patchno);
        fc2d_clawpack46_set_capacity(glob,patch,blockno,patchno);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_setaux(glob,patch,blockno,patchno);
        fc2d_clawpack5_set_capacity(glob,patch,blockno,patchno);
    }
}


void torus_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &torus_problem_setup;  /* Version-independent */

    /* Torus uses special patch setup (not advection_patch_setup), since the
       streamfunction depends on computational coordinates, not physical
       coordinates */
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    patch_vt->setup   = &torus_patch_setup;

    const user_options_t *user =  torus_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

        claw46_vt->fort_qinit  = &CLAWPACK46_QINIT;
        claw46_vt->fort_setaux = &TORUS46_SETAUX;
        claw46_vt->fort_rpn2   = &CLAWPACK46_RPN2ADV_MANIFOLD;
        claw46_vt->fort_rpt2   = &CLAWPACK46_RPT2ADV_MANIFOLD;

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt();
        
        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        claw5_vt->fort_setaux    = &TORUS5_SETAUX;
        claw5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        claw5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;
    }
}

