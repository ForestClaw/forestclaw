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

#include "filament_user.h"


void filament_problem_setup(fclaw_global_t *glob)
{
    const user_options_t* user = filament_get_options(glob);
    fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",fclaw_opt->manifold,"\% manifold\n");
        fprintf(f,  "%-24.4f   %s",user->alpha,"\% alpha\n");
        fprintf(f,  "%-24.4f   %s",user->center[0],"\% center_x\n");
        fprintf(f,  "%-24.4f   %s",user->center[1],"\% center_y\n");
        fprintf(f,  "%-24.4f   %s",user->maxelev,"\% max_elevation\n");
        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    SETPROB();
}

static
void filament_patch_setup(fclaw_global_t *glob,
                          fclaw_patch_t *patch,
                          int blockno,
                          int patchno)
{
    const user_options_t* user = filament_get_options(glob);
    claw3_advection_patch_setup_manifold(glob,patch,blockno,patchno,
                                   user->claw_version); 
}


void filament_link_solvers(fclaw_global_t *glob)
{
    fclaw_vtable_t *fc_vt = fclaw_vt(glob);
    fc_vt->problem_setup = filament_problem_setup;

    fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
    const user_options_t* user = filament_get_options(glob);

    /* example of how to set up a user defined criteria */
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->d3->fort_user_exceeds_threshold = &FCLAW3DX_USER_EXCEEDS_TH;
    
    if (user->claw_version == 4)
    {        
        if (!user->use_claw3d)
        {
            fclaw_global_essentialf("filament_user.cpp: " \
                                    "2d example not implemented.\n");
            exit(0);
        }
        else
        {
            fc3d_clawpack46_vtable_t *clawpack46_vt = fc3d_clawpack46_vt(glob);
            fc3d_clawpack46_options_t *clawopt = fc3d_clawpack46_get_options(glob);

            if (fclaw_opt->manifold)
            {
                /* This calls a manifold version of setaux */
                FCLAW_ASSERT(clawopt->mcapa != 0);
                fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
                patch_vt->setup = filament_patch_setup;
            }
            else
            {
                /* This calls a setaux routine */
                FCLAW_ASSERT(clawopt->mcapa == 0);
                clawpack46_vt->fort_setaux     = &CLAWPACK46_SETAUX;  
            }
            clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;

            /* These will work for both non-manifold and manifold case */
            clawpack46_vt->fort_rpn3       = &CLAWPACK46_RPN3;
            clawpack46_vt->fort_rpt3       = &CLAWPACK46_RPT3;
            clawpack46_vt->fort_rptt3      = &CLAWPACK46_RPTT3;            
        }
    }
    else if (user->claw_version == 5)
    {
        fclaw_global_essentialf("filament_user.cpp : " \
                                "Example not implemented for Claw version 5.\n");
        exit(0);
    }
}