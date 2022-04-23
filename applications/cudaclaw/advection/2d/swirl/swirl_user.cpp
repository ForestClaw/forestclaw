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
#include <fc2d_clawpack46.h>  
#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack46_fort.h>  
#include <clawpack46_user_fort.h>  
#include <fclaw2d_clawpatch46_fort.h>
#include "../../../../clawpack/advection/2d/all/advection_user_fort.h"

static
void swirl_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = swirl_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->period,"\% period");
        fclose(f);
    }

    /* Make sure node 0 has written 'setprob.data' before proceeding */
    fclaw2d_domain_barrier (glob->domain);
    if (user->cuda == 1)
    {
        setprob_cuda();
    }
    else{
        SETPROB();
    }
}


void swirl_link_solvers(fclaw2d_global_t *glob)
{
	fclaw2d_vtable_t *vt = fclaw2d_vt();
    //vt->problem_setup = &swirl_problem_setup;  /* Version-independent */
    //fclaw2d_patch_vtable_t*  patch_vt = fclaw2d_patch_vt();  
    const user_options_t* user = swirl_get_options(glob);

    if (user->cuda == 0)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();        
        
        clawpack46_vt->fort_qinit     = &CUDACLAW_QINIT;
        clawpack46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
        clawpack46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;

        /* Velocity is set here rather than in setaux, because we have a 
           time dependent velocity field */
        clawpack46_vt->fort_b4step2   = &CLAWPACK46_B4STEP2;
    }
    else
    {
        vt->problem_setup = &swirl_problem_setup;  /* Version-independent */

        // const user_options_t* user = swirl_get_options(glob);
        fc2d_cudaclaw_vtable_t *cudaclaw_vt = fc2d_cudaclaw_vt();        

        cudaclaw_vt->fort_qinit     = &CUDACLAW_QINIT;
            
        //cudaclaw_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
        swirl_assign_rpn2(&cudaclaw_vt->cuda_rpn2);
        FCLAW_ASSERT(cudaclaw_vt->cuda_rpn2 != NULL);

        // cudaclaw_vt->fort_b4step2   = &CUDACLAW_B4STEP2;
        swirl_assign_b4step2(&cudaclaw_vt->cuda_b4step2);
        FCLAW_ASSERT(cudaclaw_vt->cuda_b4step2 != NULL);

        //cudaclaw_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;
        swirl_assign_rpt2(&cudaclaw_vt->cuda_rpt2);
        FCLAW_ASSERT(cudaclaw_vt->cuda_rpt2 != NULL);
    }
}




