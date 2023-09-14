/*
Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright

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

#include "phasefield_user.h"
#include "phasefield_options.h"

#include "phasefield_operator.h"

#include <fclaw_include_all.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_fort.h>
#include <fc2d_thunderegg_options.h>
#include <fc2d_thunderegg_physical_bc.h>
#include <fc2d_thunderegg_starpatch.h>
#include <fc2d_thunderegg_fivepoint.h>

#include <fclaw_elliptic_solver.h>


#include <fclaw_farraybox.hpp>


static
void phasefield_problem_setup(fclaw_global_t *glob)
{
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        const phasefield_options_t* user = phasefield_get_options(glob);

        fprintf(f,  "%-24d   %s",  user->example,    "% example\n");
        fprintf(f,  "%-24.6f   %s",user->S,          "% S\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,      "% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->m,          "% m\n");
        fprintf(f,  "%-24.6f   %s",user->xi,         "% xi\n");
        fprintf(f,  "%-24.6f   %s",user->k,          "% k\n");
        fprintf(f,  "%-24.6f   %s",user->gamma,      "% gamma\n");
        fprintf(f,  "%-24.6f   %s",user->r0,         "% r0\n");
        fprintf(f,  "%-24.6f   %s",user->x0,         "% x0\n");
        fprintf(f,  "%-24.6f   %s",user->y0,         "% y0\n");

#if 0
        const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
        double xlower = fclaw_opt->ax;
        double xupper = fclaw_opt->bx;
        double ylower = fclaw_opt->ay;
        double yupper = fclaw_opt->by;
#endif        

        /* These are passed to fortran boundary condition routines */
        fc2d_thunderegg_options_t*  mg_opt = fc2d_thunderegg_get_options(glob);    
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[0],  "% bc[0]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[1],  "% bc[1]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[2],  "% bc[2]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[3],  "% bc[3]\n");


        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    PHASEFIELD_SETPROB(); /* This file reads the file just created above */
}

static
void phasefield_initialize(fclaw_global_t *glob,
                           fclaw_patch_t *patch,
                           int blockno,
                           int patchno)
{

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    PHASEFIELD_INIT(&meqn, &mbc, &mx, &my, &xlower, &ylower, &dx, &dy,q);
}


static
void phasefield_rhs(fclaw_global_t *glob,
                fclaw_patch_t *patch,
                int blockno,
                int patchno)
{

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
    fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    /* This function supplies an analytic right hand side. */
    int method = 1;
    double dt = glob->curr_dt;
    PHASEFIELD_FORT_RHS(&blockno, &mbc, &mx, &my, &meqn, &mfields,
                  &xlower, &ylower, &dx, &dy,&dt, &method,q,rhs);
}



#if 0
static
void phasefield_time_header_ascii(fclaw_global_t* glob, int iframe)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                fclaw_clawpatch_get_options(glob);
    char matname1[20];
    sprintf(matname1,"fort.q%04d",iframe);

#if 1
    FILE *f1 = fopen(matname1,"w");
    fclose(f1);
#endif    


    char matname2[20];
    sprintf(matname2,"fort.t%04d",iframe);

    double time = glob->curr_time;

    int ngrids = glob->domain->global_num_patches;

    int mfields = clawpatch_opt->rhs_fields;  
    int maux = clawpatch_opt->maux;


    FILE *f2 = fopen(matname2,"w");
    fprintf(f2,"%12.6f %23s\n%5d %30s\n%5d %30s\n%5d %30s\n%5d %30s\n",time,"time",
            mfields+2,"mfields",ngrids,"ngrids",maux,"num_aux",2,"num_dim");
    fclose(f2);

#if 0
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    /* header writes out mfields+2 fields (computed soln, true soln, error); */
    clawpatch_vt->fort_header_ascii(matname1,matname2,&time,&mfields,&maux,&ngrids);
#endif    

}


static
void cb_phasefield_output_ascii(fclaw_domain_t * domain,
                            fclaw_patch_t * patch,
                            int blockno, int patchno,
                            void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_global_t *glob = (fclaw_global_t*) s->glob;
    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    int global_num, local_num;
    int level;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,&local_num, &level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    double *error = fclaw2d_clawpatch_get_error(glob,patch);
    double *soln  = fclaw2d_clawpatch_get_exactsoln(glob,patch);

    char fname[BUFSIZ];
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);

    /* The fort routine is defined by a clawpack solver and handles 
       the layout of q in memory (i,j,m) or (m,i,j), etc */
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_output_ascii);

    PHASEFIELD_FORT_OUTPUT_ASCII(fname,&mx,&my,&mfields,&mbc,
                             &xlower,&ylower,&dx,&dy,rhs,
                             soln, error,
                             &global_num,&level,&blockno,
                             &glob->mpirank);
}
#endif


int phasefield_tag4refinement(fclaw_global_t *glob,
                             fclaw_patch_t *this_patch,
                             int blockno, int patchno,
                             int initflag)
{
    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    int tag_patch;
    double refine_threshold;

    refine_threshold = fclaw_opt->refine_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    tag_patch = 0;
    clawpatch_vt->d2->fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                          &blockno, q,&refine_threshold,
                                          &initflag,&tag_patch);
    return tag_patch;
}

static
int phasefield_tag4coarsening(fclaw_global_t *glob,
                             fclaw_patch_t *fine_patches,
                             int blockno,
                             int patchno,
                             int initflag)
{
    fclaw_patch_t *patch0 = &fine_patches[0];

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q[4];
    int meqn;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
    }

    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

    int tag_patch = 0;
    clawpatch_vt->d2->fort_tag4coarsening(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                          &blockno, q[0],q[1],q[2],q[3],
                                          &coarsen_threshold,&initflag,&tag_patch);
    return tag_patch == 1;
}

static
void phasefield_bc2(fclaw_global_t *glob,
                   fclaw_patch_t *patch,
                   int block_idx,
                   int patch_idx,
                   double t,
                   double dt,
                   int intersects_bc[],
                   int time_interp)
{


    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    PHASEFIELD_FORT_BC2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                           &dx,&dy,q,&t,&dt,intersects_bc);
}


void phasefield_link_solvers(fclaw_global_t *glob)
{
#if 0 
    /* These are listed here for reference */
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt(glob);
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt(glob);
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
#endif
    /* ForestClaw vtable */
    fclaw_vtable_t *fc_vt = fclaw_vt(glob);
    fc_vt->problem_setup = &phasefield_problem_setup;  

    /* Patch : RHS function */
    fclaw_patch_vtable_t* patch_vt = fclaw_patch_vt(glob);
    patch_vt->physical_bc = phasefield_bc2;   /* Doesn't do anything */
    patch_vt->rhs = phasefield_rhs;          /* Overwrites default */
    patch_vt->initialize = phasefield_initialize;   /* Get an initial refinement */

    /* Assign phasefield operator vtable */
    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);
    mg_vt->patch_operator = phasefield_solve;

    mg_vt->fort_apply_bc = &PHASEFIELD_FORT_APPLY_BC;
    mg_vt->fort_eval_bc  = &PHASEFIELD_NEUMANN;   // For non-homogeneous BCs

    // tagging routines
    patch_vt->tag4refinement       = phasefield_tag4refinement;
    patch_vt->tag4coarsening       = phasefield_tag4coarsening;

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->d2->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->d2->fort_tag4coarsening = &TAG4COARSENING;

}

