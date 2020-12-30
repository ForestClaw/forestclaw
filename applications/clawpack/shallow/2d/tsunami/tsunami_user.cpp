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

#include "tsunami_user.h"

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>

#include <clawpack46_user_fort.h>

void tsunami_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    vt->problem_setup = &tsunami_problem_setup;  /* Version-independent */

    //fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

    claw46_vt->fort_qinit  = &CLAWPACK46_QINIT;
    claw46_vt->fort_setaux = &CLAWPACK46_SETAUX;
    claw46_vt->fort_rpn2   = &RPN2_TSUNAMI;
    claw46_vt->fort_rpt2   = &RPT2_TSUNAMI;
}


void tsunami_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = tsunami_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->g,"\% grav");
        fprintf(f,"%-24.4f %s\n",user->a,"\% a");
        fprintf(f,"%-24.4f %s\n",user->b,"\% b");
        fprintf(f,"%-24.4f %s\n",user->h0,"\% h0");
        fprintf(f,"%-24.4f %s\n",user->breaking,"\% breaking");
        fprintf(f,"%-24.4f %s\n",user->alpha,"\% alpha");
        fprintf(f,"%-24.4f %s\n",user->dry_tolerance,"\% dry_tolerance");
        fprintf(f,"%-24.4f %s\n",user->sea_level,"\% sea_level");        
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    TSUNAMI_SETPROB();  /* Reads file created above */
}


