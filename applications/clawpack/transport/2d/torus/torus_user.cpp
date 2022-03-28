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
    transport_patch_setup_manifold(glob,patch,blockno,patchno,
                                   user->claw_version);

    /* Torus velocity field is not time dependent, so we can set up the 
       velocity here, using b4step */

    double t = 0;
    double dt = -1;    /* Not used */
    transport_b4step2_manifold(glob,patch,blockno,patchno,t, dt,
                               user->claw_version);
}


void torus_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &torus_problem_setup;  /* Version-independent */

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
    patch_vt->setup   = &torus_patch_setup;

    const user_options_t *user =  torus_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

        claw46_vt->fort_qinit = &CLAWPACK46_QINIT;

        fc2d_clawpack46_options_t *claw46_opt = fc2d_clawpack46_get_options(glob);

        /* Solve conservative equation using cell-centered velocity */
        /* Fwaves give best accuracy */
        claw46_opt->use_fwaves = 1;
        claw46_vt->fort_rpn2fw = CLAWPACK46_RPN2CONS_FW_MANIFOLD; 

        /* Transverse solver : Conservative form */
        claw46_vt->fort_rpt2fw  = &CLAWPACK46_RPT2CONS_MANIFOLD;  

        /* Flux function (for conservative fix) */
         claw46_vt->fort_rpn2_cons = &RPN2CONS_UPDATE_MANIFOLD;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt();

        claw5_vt->fort_qinit = &CLAWPACK5_QINIT;

        fc2d_clawpack5_options_t *claw5_opt = fc2d_clawpack5_get_options(glob);

        /* Solve conservative equation using cell-centered velocity */
        /* Fwaves give best accuracy */
        claw5_opt->use_fwaves = 1;
        claw5_vt->fort_rpn2 = CLAWPACK5_RPN2CONS_MANIFOLD; 

        /* Transverse solver : Conservative form */
        claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2CONS_MANIFOLD;  

        /* Flux function (for conservative fix) */
         claw5_vt->fort_rpn2_cons = &RPN2CONS_UPDATE_MANIFOLD;
    }
}

