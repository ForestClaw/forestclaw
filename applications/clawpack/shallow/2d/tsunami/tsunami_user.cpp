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

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fc2d_clawpack46.h>
#include <clawpack46_user_fort.h>

#include <fclaw2d_elliptic_solver.h>

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_fort.h>
#include <fc2d_thunderegg_options.h>
#include <fc2d_thunderegg_physical_bc.h>
#include <fc2d_thunderegg_starpatch.h>
#include <fc2d_thunderegg_fivepoint.h>

#include "../sgn/sgn_options.h"
#include "../sgn/sgn_operator.h"
#include "../sgn/sgn_patch_operator.h"

#include "../sgn/sgn_fort.h"



static
void tsunami_problem_setup(fclaw2d_global_t* glob)
{
    const sgn_options_t* sgn = sgn_get_options(glob);
    const user_options_t* user = tsunami_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->g,"\% grav");
        fprintf(f,"%-24.4f %s\n",user->a,"\% a");
        fprintf(f,"%-24.4f %s\n",user->b,"\% b");
        fprintf(f,"%-24.4f %s\n",user->h0,"\% h0");

        fprintf(f,"%-24.4f %s\n",sgn->breaking,"\% breaking");
        fprintf(f,"%-24.4f %s\n",sgn->alpha,"\% alpha");
        fprintf(f,"%-24.4f %s\n",sgn->dry_tolerance,"\% dry_tolerance");
        fprintf(f,"%-24.4f %s\n",sgn->sea_level,"\% sea_level");        
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    TSUNAMI_SETPROB();  /* Reads file created above */
}


static
void tsunami_rhs(fclaw2d_global_t *glob,
                 fclaw2d_patch_t *patch,
                 int blockno,
                 int patchno)
{

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    int meqn;
    double *q;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    SGN_FORT_RHS(&blockno, &mbc, &mx, &my, &meqn, &mfields,
                 &xlower, &ylower, &dx, &dy,q,rhs);
}

#if 0
int tsunami_tag4refinement(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *this_patch,
                             int blockno, int patchno,
                             int initflag)
{
    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    int tag_patch;
    double refine_threshold;

    refine_threshold = fclaw_opt->refine_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    tag_patch = 0;
    clawpatch_vt->fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                      &blockno, q,&refine_threshold,
                                      &initflag,&tag_patch);
    return tag_patch;
}

static
int tsunami_tag4coarsening(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *fine_patches,
                             int blockno,
                             int patchno,
                             int initflag)
{
    fclaw2d_patch_t *patch0 = &fine_patches[0];

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q[4];
    int meqn;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
    }

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();

    int tag_patch = 0;
    clawpatch_vt->fort_tag4coarsening(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                      &blockno, q[0],q[1],q[2],q[3],
                                      &coarsen_threshold,&initflag,&tag_patch);
    return tag_patch == 1;
}
#endif

void tsunami_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    vt->problem_setup = &tsunami_problem_setup;  /* Version-independent */

    //fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

    claw46_vt->fort_qinit  = &CLAWPACK46_QINIT;
    claw46_vt->fort_setaux = &CLAWPACK46_SETAUX;
    claw46_vt->fort_rpn2   = &RPN2_TSUNAMI;     /* or RPN2_GEOCLAW; */
    claw46_vt->fort_rpt2   = &RPT2_TSUNAMI;     /* or RPT2_GEOCLAW; */

    /* for SGN operator */
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    patch_vt->rhs = tsunami_rhs;          /* Overwrites default */

    /* Assign SGN operator vtable */
    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt();
    mg_vt->patch_operator = sgn_solve;

    mg_vt->fort_apply_bc = &SGN_FORT_APPLY_BC;
    mg_vt->fort_eval_bc  = &SGN_NEUMANN;   // For non-homogeneous BCs

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;
}





