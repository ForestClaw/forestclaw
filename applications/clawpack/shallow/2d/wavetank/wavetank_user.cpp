/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "wavetank_user.h"
#include "../sgn/sgn.h"
#include "../sgn/sgn_fort.h"

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


#include <clawpack46_user_fort.h>

static
void wavetank_problem_setup(fclaw2d_global_t* glob)
{
    const sgn_options_t* sgn = sgn_get_options(glob);
    const user_options_t* user = wavetank_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->g,"\% grav");
        fprintf(f,"%-24.4f %s\n",user->dry_tolerance,"\% dry_tolerance");
        fprintf(f,"%-24.4f %s\n",user->sea_level,"\% sea_level");        

        fprintf(f,"%-24.4f %s\n",sgn->breaking,"\% breaking");
        fprintf(f,"%-24.4f %s\n",sgn->alpha,"\% alpha");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    WAVETANK_SETPROB();  /* Reads file created above */
}


/* --------------------------------- Output functions ---------------------------- */

#if 1
int wavetank_tag4refinement(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *patch,
                             int blockno, int patchno,
                             int initflag)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    int tag_patch;
    double refine_threshold;

    refine_threshold = fclaw_opt->refine_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    tag_patch = 0;
    int level = patch->level;
    double time = glob->curr_time;
    WAVETANK_FORT_TAG4REFINEMENT(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                 &time, &blockno, q,&refine_threshold,
                                 &level,
                                 &initflag,&tag_patch);
    return tag_patch;
}
#endif

#if 1
static
int wavetank_tag4coarsening(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *fine_patches,
                             int blockno,
                             int patchno,
                             int initflag)
{
    fclaw2d_patch_t *patch0 = &fine_patches[0];

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int mx,my,mbc;
    double xlower[4],ylower[4],dx,dy;
    for(int k = 0; k < 4; k++)
    {
        fclaw2d_clawpatch_grid_data(glob,patch0+k,&mx,&my,&mbc,
                                    &xlower[k],&ylower[k],&dx,&dy);        
    }

    double *q[4];
    int meqn;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
    }

    int tag_patch = 0;
    double time = glob->curr_time;
    WAVETANK_FORT_TAG4COARSENING(&mx,&my,&mbc,&meqn,xlower,ylower,&dx,&dy,
                                 &time, 
                                 &blockno, q[0],q[1],q[2],q[3],
                                 &coarsen_threshold,&initflag,&tag_patch);

#if 1
    if (time > 15 && tag_patch == 1) 
    {
        printf("coarsen %d\n",tag_patch);

    }
#endif    
    return tag_patch == 1;
}
#endif

void wavetank_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();

    fclaw_vt->problem_setup = &wavetank_problem_setup;  /* Version-independent */
    fclaw_vt->output_frame  = fc2d_geoclaw_output_ascii;

    //fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

    claw46_vt->fort_qinit  = &CLAWPACK46_QINIT;
    claw46_vt->fort_setaux = &CLAWPACK46_SETAUX;
    claw46_vt->fort_rpn2   = &RPN2_GEOCLAW;     /* or RPN2_GEOCLAW; */
    claw46_vt->fort_rpt2   = &RPT2_GEOCLAW;     /* or RPT2_GEOCLAW; */
    claw46_vt->fort_bc2    = &CLAWPACK46_BC2;     /* or RPT2_GEOCLAW; */


    /* for SGN operator */
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    patch_vt->rhs = sgn_rhs;          /* Overwrites default */

    /* Assign SGN operator vtable */
    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt();
    mg_vt->patch_operator = sgn_solve;

    mg_vt->fort_apply_bc = &SGN_FORT_APPLY_BC;
    mg_vt->fort_eval_bc  = &SGN_NEUMANN;   // For non-homogeneous BCs

    patch_vt->tag4refinement = wavetank_tag4refinement;
    patch_vt->tag4coarsening = wavetank_tag4coarsening;
#if 0
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;
#endif    

}





