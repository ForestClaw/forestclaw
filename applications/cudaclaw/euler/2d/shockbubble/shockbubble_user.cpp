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

#include "shockbubble_user.h"
#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>  

#include "../../../../clawpack/euler/2d/rp/euler_user_fort.h"

void shockbubble_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = shockbubble_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24.16f   %s",user->gamma,"\% gamma\n");
        fprintf(f,  "%-24.16f   %s",user->x0,"\% x0\n");
        fprintf(f,  "%-24.16f   %s",user->y0,"\% y0\n");
        fprintf(f,  "%-24.16f   %s",user->r0,"\% r0\n");
        fprintf(f,  "%-24.16f   %s",user->rhoin,"\% rhoin\n");
        fprintf(f,  "%-24.16f   %s",user->pinf,"\% pinf\n");
        fprintf(f,  "%-24d   %s",   user->idisc,"\% idisc\n");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
    fclaw2d_domain_barrier (glob->domain);  
    if (user->cuda != 0)
    {
        /* Call CUDA setprob to set parameters needed for Riemann solvers */ 
        setprob_cuda(); 
    }
    /* Call Fortran setprob to get parameters needed for qinit */
    SETPROB();
 
    
}

#if 0
static
void shockbubble_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = shockbubble_get_options(glob);

    if (glob->mpirank == 0)
    {

    }

    SHOCKBUBBLE_SETPROB(&user->gamma, &user->x0, &user->y0, &user->r0,
                        &user->rhoin, &user->pinf, &user->idisc);
    if (user->cuda == 1)
    {
        shockbubble_setprob_cuda(user->gamma);
    }

}
#endif

void shockbubble_link_solvers(fclaw2d_global_t *glob)
{
    // const user_options_t* user = shockbubble_get_options(glob);
    // fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

    // fc2d_clawpack46_options_t *clawopt = fc2d_clawpack46_get_options(glob);

    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    const user_options_t* user = shockbubble_get_options(glob);
    fclaw_vt->problem_setup = &shockbubble_problem_setup;


    if (user->cuda == 0)
    {
        
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();
        fc2d_clawpack46_options_t *clawopt = fc2d_clawpack46_get_options(glob);

        claw46_vt->fort_qinit  = &CUDACLAW_QINIT;
        claw46_vt->fort_setaux = &CLAWPACK46_SETAUX;
        claw46_vt->fort_bc2    = &CUDACLAW_BC2;   /* Special  BCs at left edge */
        claw46_vt->fort_src2   = &CUDACLAW_SRC2;  /* To simulate axis-symmetric */

        // switch (clawopt->mwaves)
        // {
        // case 4:
        //     /* Requires meqn=4 */
            claw46_vt->fort_rpn2   = &CLAWPACK46_RPN2_EULER4;  /* No tracer */
            claw46_vt->fort_rpt2   = &CLAWPACK46_RPT2_EULER4;
        //     break;
        // case 5:
        //     /* Requires meqn=5 */
        //     claw46_vt->fort_rpn2   = &CLAWPACK46_RPN2_EULER5;  /* Includes a tracer */
        //     claw46_vt->fort_rpt2   = &CLAWPACK46_RPT2_EULER5;
        //     break;
        // default:
            // SC_ABORT_NOT_REACHED ();
        // }
    }
    else
    {

        //fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();
        fc2d_cudaclaw_vtable_t *cuclaw_vt   = fc2d_cudaclaw_vt();
        cuclaw_vt->fort_qinit  = &CUDACLAW_QINIT;
        cuclaw_vt->fort_bc2    = &CUDACLAW_BC2;   /* Special  BCs at left edge */
        cuclaw_vt->fort_setaux = &CLAWPACK46_SETAUX;
        cuclaw_vt->fort_src2   = &CUDACLAW_SRC2;  /* To simulate axis-symmetric */

        shockbubble_assign_rpn2(&cuclaw_vt->cuda_rpn2);
        FCLAW_ASSERT(cuclaw_vt->cuda_rpn2 != NULL);

        shockbubble_assign_rpt2(&cuclaw_vt->cuda_rpt2);
        FCLAW_ASSERT(cuclaw_vt->cuda_rpt2 != NULL);
    }
#if 0
    /* Use divided differences to tag grids */
    clawpatch_vt->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;
#endif    
}