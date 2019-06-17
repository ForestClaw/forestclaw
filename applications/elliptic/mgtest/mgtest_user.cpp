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

#include <fclaw2d_elliptic_solver.h>


static
void mgtest_rhs(fclaw2d_global_t *glob,
                fclaw2d_patch_t *patch,
                int blockno,
                int patchno);
    

void mgtest_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    fclaw_vt->problem_setup = &mgtest_problem_setup;  
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

    /* RHS function */
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    //patch_vt->rhs = mgtest_rhs;
    patch_vt->initialize = mgtest_rhs;   /* To get an initial refinement */

    /* Only needed if Fortran subroutine is useful and can be customized */
    fc2d_multigrid_vtable_t*  mg_vt = fc2d_multigrid_vt();
    mg_vt->fort_rhs = &MGTEST_FORT_RHS;

    /* Compute the error */
    clawpatch_vt->fort_compute_patch_error = &MGTEST_COMPUTE_ERROR;


    /* Clawpatch tagging routines */
    clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;


    /* Do something with user options? */
    //const user_options_t* user_opt = mgtest_get_options(glob);
}

void mgtest_problem_setup(fclaw2d_global_t* glob)
{
    const mgtest_options_t* user = mgtest_get_options(glob);

    MGTEST_SETPROB(&user->rhs_choice, &user->alpha,
                   &user->x0, &user->y0, 
                   &user->a, &user->b);
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





