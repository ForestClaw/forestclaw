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

#include "radial_user.h"

static
void radial_problem_setup(fclaw2d_global_t* glob)
{
    user_options_t* user = radial_get_options(glob);    
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24.16f   %s",user->rho,"\% rho\n");
        fprintf(f,  "%-24.16f   %s",user->bulk,"\% bulk\n");
        fclose(f);
    }
    /* We want to make sure node 0 gets here before proceeding */
    fclaw2d_domain_barrier (glob->domain);  /* redundant?  */

    if(user->cuda != 0)
    {
        setprob_cuda();

    //SETPROB();
    }

    SETPROB();
      
}

void radial_link_solvers(fclaw2d_global_t *glob)
{
    //fclaw2d_patch_vtable_t*  patch_vt = fclaw2d_patch_vt(glob);  

    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &radial_problem_setup;  /* Version-independent */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob); //added
    
    const user_options_t* user = radial_get_options(glob);
     if(user->cuda != 0)
     {
        //const user_options_t* user = radial_get_options(glob);
        fc2d_cudaclaw_vtable_t *cudaclaw_vt = fc2d_cudaclaw_vt(glob);        
        cudaclaw_vt->fort_qinit     = &CUDACLAW_QINIT;

        radial_assign_rpn2(&cudaclaw_vt->cuda_rpn2);
        FCLAW_ASSERT(cudaclaw_vt->cuda_rpn2 != NULL);

        radial_assign_rpt2(&cudaclaw_vt->cuda_rpt2);
        FCLAW_ASSERT(cudaclaw_vt->cuda_rpt2 != NULL);
    }
    else
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);

        claw46_vt->fort_qinit     = &CUDACLAW_QINIT;
        
        claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2;
        claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2;
  
        
    }
    
}

