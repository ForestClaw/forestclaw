/*
Copyright (c) 2019-2020 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright

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

#include "allencahn_user.h"
#include "allencahn_options.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fc2d_multigrid.h>
#include <fc2d_multigrid_fort.h>
#include <fc2d_multigrid_options.h>
#include <fc2d_multigrid_physical_bc.h>
#include <fc2d_multigrid_starpatch.h>
#include <fc2d_multigrid_fivepoint.h>

#include <fclaw2d_elliptic_solver.h>


#include <fclaw2d_farraybox.hpp>


static
void allencahn_problem_setup(fclaw2d_global_t *glob)
{
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        const allencahn_options_t* user = allencahn_get_options(glob);

        fprintf(f,  "%-24d   %s",  user->example,    "% example\n");
        fprintf(f,  "%-24d   %s",  user->beta_choice,"% beta_choice\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,      "% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->x0,         "% x0\n");
        fprintf(f,  "%-24.6f   %s",user->y0,         "% y0\n");
        fprintf(f,  "%-24.6f   %s",user->a,          "% a\n");
        fprintf(f,  "%-24.6f   %s",user->b,          "% b\n");
        fprintf(f,  "%-24.6f   %s",user->eps_disk,   "% eps_disk\n");
        fprintf(f,  "%-24d   %s",user->m_polar,    "% m_polar\n");
        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% x0[%d]\n",user->x0_polar[i],i); 

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% y0[%d]\n",user->y0_polar[i],i);            

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% r0[%d]\n",user->r0_polar[i],i);            

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% r1[%d]\n",user->r1_polar[i],i);            

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24d   %% n[%d]\n",user->n_polar[i],i);            

        fc2d_multigrid_options_t*  mg_opt = fc2d_multigrid_get_options(glob);    
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[0],  "% bc[0]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[1],  "% bc[1]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[2],  "% bc[2]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[3],  "% bc[3]\n");


        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    ALLENCAHN_SETPROB(); /* This file reads the file just created above */
}

static
void allencahn_initialize(fclaw2d_global_t *glob,
                     fclaw2d_patch_t *patch,
                     int blockno,
                     int patchno)
{

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *q;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    /* This function supplies an analytic right hand side. */

    ALLENCAHN_INIT(&meqn, &mbc, &mx, &my, &xlower, &ylower, &dx, &dy,q);

}


static
void allencahn_rhs(fclaw2d_global_t *glob,
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

    /* This function supplies an analytic right hand side. */
    int method = 1;
    double dt = glob->curr_dt;
    ALLENCAHN_FORT_RHS(&blockno, &mbc, &mx, &my, &meqn, &mfields,
                  &xlower, &ylower, &dx, &dy,&dt, &method,q,rhs);

}



#if 0
static
void allencahn_time_header_ascii(fclaw2d_global_t* glob, int iframe)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
                fclaw2d_clawpatch_get_options(glob);
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
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    /* header writes out mfields+2 fields (computed soln, true soln, error); */
    clawpatch_vt->fort_header_ascii(matname1,matname2,&time,&mfields,&maux,&ngrids);
#endif    

}


static
void cb_allencahn_output_ascii(fclaw2d_domain_t * domain,
                            fclaw2d_patch_t * patch,
                            int blockno, int patchno,
                            void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;
    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    int global_num, local_num;
    int level;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,&local_num, &level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    double *error = fclaw2d_clawpatch_get_error(glob,patch);
    double *soln  = fclaw2d_clawpatch_get_exactsoln(glob,patch);

    char fname[BUFSIZ];
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);

    /* The fort routine is defined by a clawpack solver and handles 
       the layout of q in memory (i,j,m) or (m,i,j), etc */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    FCLAW_ASSERT(clawpatch_vt->fort_output_ascii);

    ALLENCAHN_FORT_OUTPUT_ASCII(fname,&mx,&my,&mfields,&mbc,
                             &xlower,&ylower,&dx,&dy,rhs,
                             soln, error,
                             &global_num,&level,&blockno,
                             &glob->mpirank);
}
#endif


int allencahn_tag4refinement(fclaw2d_global_t *glob,
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

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,this_patch,&rhs,&mfields);

    tag_patch = 0;
    clawpatch_vt->fort_tag4refinement(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs,&refine_threshold,
                                      &initflag,&tag_patch);
    return tag_patch;
}

static
int allencahn_tag4coarsening(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *fine_patches,
                             int blockno,
                             int patchno)
{
    fclaw2d_patch_t *patch0 = &fine_patches[0];

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs[4];
    int mfields;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_rhs_data(glob,&fine_patches[igrid],&rhs[igrid],&mfields);
    }

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();

    int tag_patch = 0;
    clawpatch_vt->fort_tag4coarsening(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs[0],rhs[1],rhs[2],rhs[3],
                                      &coarsen_threshold,&tag_patch);
    return tag_patch == 1;
}


void allencahn_link_solvers(fclaw2d_global_t *glob)
{
#if 0 
    /* These are listed here for reference */
    fc2d_multigrid_options_t *mg_opt = fc2d_multigrid_get_options(glob);
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
#endif
    /* ForestClaw vtable */
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    fclaw_vt->problem_setup = &allencahn_problem_setup;  

    /* Patch : RHS function */
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    patch_vt->physical_bc = fclaw2d_physical_bc_default;   /* Doesn't do anything */
    patch_vt->rhs = allencahn_rhs;          /* Overwrites default */
    patch_vt->initialize = allencahn_initialize;   /* Get an initial refinement */

    /* Multigrid vtable */
    fc2d_multigrid_vtable_t*  mg_vt = fc2d_multigrid_vt();
    //mg_vt->fort_rhs       = &ALLENCAHN_FORT_RHS;

    mg_vt->fort_beta      = &ALLENCAHN_FORT_BETA;
    
    mg_vt->fort_apply_bc = &ALLENCAHN_FORT_APPLY_BC;
    mg_vt->fort_eval_bc  = &ALLENCAHN_DIRICHLET;   // For non-homogeneous BCs

    // tagging routines
    patch_vt->tag4refinement       = allencahn_tag4refinement;
    patch_vt->tag4coarsening       = allencahn_tag4coarsening;

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;

}

