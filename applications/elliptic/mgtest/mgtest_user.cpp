/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright

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

#include "mgtest_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fc2d_multigrid.h>
#include <fc2d_multigrid_fort.h>
#include <fc2d_multigrid_options.h>
#include <fc2d_multigrid_physical_bc.h>

#include <fclaw2d_elliptic_solver.h>


static
void mgtest_problem_setup(fclaw2d_global_t *glob)
{
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        const mgtest_options_t* user = mgtest_get_options(glob);

        fprintf(f,  "%-24d   %s",  user->example,"\% example\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,  "\% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->x0,     "\% x0\n");
        fprintf(f,  "%-24.6f   %s",user->y0,     "\% y0\n");
        fprintf(f,  "%-24.6f   %s",user->a,      "\% a\n");
        fprintf(f,  "%-24.6f   %s",user->b,      "\% b\n");

        fc2d_multigrid_options_t*  mg_opt = fc2d_multigrid_get_options(glob);    
        fprintf(f,  "%24d      %s",mg_opt->boundary_conditions[0],  "\% bc[0]\n");
        fprintf(f,  "%24d      %s",mg_opt->boundary_conditions[1],  "\% bc[1]\n");
        fprintf(f,  "%24d      %s",mg_opt->boundary_conditions[2],  "\% bc[2]\n");
        fprintf(f,  "%24d      %s",mg_opt->boundary_conditions[3],  "\% bc[3]\n");

        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    MGTEST_SETPROB(); /* This file reads the file just created above */
}


static
void mgtest_rhs(fclaw2d_global_t *glob,
                fclaw2d_patch_t *patch,
                int blockno,
                int patchno)
{
    int mx,my,mbc,meqn;
    double dx,dy,xlower,ylower;
    double *q;

    fc2d_multigrid_vtable_t*  mg_vt = fc2d_multigrid_vt();

    /* Compute right hand side */
    FCLAW_ASSERT(mg_vt->fort_rhs != NULL);

    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    /* Or some other suitable function that sets up rhs on patches */
    mg_vt->fort_rhs(&blockno,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q);
}


void mgtest_link_solvers(fclaw2d_global_t *glob)
{
    /* ForestClaw vtable */
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    fclaw_vt->problem_setup = &mgtest_problem_setup;  

    /* Patch : RHS function */
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    patch_vt->rhs = mgtest_rhs;          /* Overwrites default */
    patch_vt->initialize = mgtest_rhs;   /* Get an initial refinement */

    /* Multigrid vtable */
    fc2d_multigrid_vtable_t*  mg_vt = fc2d_multigrid_vt();
    mg_vt->fort_rhs      = &MGTEST_FORT_RHS;
    mg_vt->fort_apply_bc = &MGTEST_FORT_APPLY_BC;
    mg_vt->fort_eval_bc  = &MGTEST_FORT_EVAL_BC;   // For non-homogeneous BCs

    /* Clawpatch : Compute the error */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    clawpatch_vt->fort_compute_patch_error = &MGTEST_COMPUTE_ERROR;

    // tagging routines
    clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;

    /* Do something with user options? */
    //const user_options_t* user_opt = mgtest_get_options(glob);
}

