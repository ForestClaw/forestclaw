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

#include "poisson_user.h"


static
void poisson_problem_setup(fclaw_global_t *glob)
{
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        const poisson_options_t* user = poisson_get_options(glob);

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

        fc2d_thunderegg_options_t*  mg_opt = fc2d_thunderegg_get_options(glob);    
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[0],  "% bc[0]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[1],  "% bc[1]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[2],  "% bc[2]\n");
        fprintf(f,  "%-24d   %s",mg_opt->boundary_conditions[3],  "% bc[3]\n");


        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    POISSON_SETPROB(); /* This file reads the file just created above */
}


static
void poisson_rhs(fclaw_global_t *glob,
                fclaw_patch_t *patch,
                int blockno,
                int patchno)
{

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
    fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    /* Compute right hand side */
    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);
    FCLAW_ASSERT(mg_vt->fort_rhs != NULL);

    /* This function supplies an analytic right hand side. */
    mg_vt->fort_rhs(&blockno,&mbc,&mx,&my,&mfields, &xlower,&ylower,&dx,&dy,rhs);
}

static void poisson_patch_setup(struct fclaw_global *glob,
                                struct fclaw_patch *patch,
                                int blockno,
                                int patchno)
{
    poisson_rhs(glob,patch,blockno,patchno);
}

static
void poisson_compute_error(fclaw_global_t *glob,
                          fclaw_patch_t *patch,
                          int blockno,
                          int patchno,
                          void *user)
{
    poisson_error_info_t* error_data = (poisson_error_info_t*) user;
    //const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    if (clawpatch_vt->d2->fort_compute_patch_error != NULL)
    {
        int mx, my, mbc;
        double xlower,ylower,dx,dy;
        fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,&xlower,
                                    &ylower,&dx,&dy);

        double *area = fclaw2d_clawpatch_get_area(glob,patch);  /* Might be null */

        /* Solution is stored in the RHS */
        double *rhs, *err, *soln;  
        int mfields;
        fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);
        fclaw_clawpatch_elliptic_error_data(glob,patch,&err,&mfields);
        fclaw_clawpatch_elliptic_soln_data(glob,patch,&soln,&mfields);

        double t = glob->curr_time;

        clawpatch_vt->d2->fort_compute_patch_error(&blockno, &mx,&my,&mbc,
                                                   &mfields,&dx,&dy,
                                                   &xlower,&ylower, &t, rhs, err, soln);

        /* Accumulate sums and maximums needed to compute error norms */

        FCLAW_ASSERT(clawpatch_vt->d2->fort_compute_error_norm != NULL);
        clawpatch_vt->d2->fort_compute_error_norm(&blockno, &mx, &my, &mbc, &mfields, 
                                                  &dx,&dy, area, err,
                                                  error_data->local_error);

    }
}


static
void poisson_conservation_check(fclaw_global_t *glob,
                               fclaw_patch_t *patch,
                               int blockno,
                               int patchno,
                               void *user)
{
    poisson_error_info_t* error_data = (poisson_error_info_t*) user;
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;  /* Solution is stored in the right hand side */ 
    fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->d2->fort_conservation_check != NULL);


    /* Need a better way to determine which diagnostic to do */
    double* area = fclaw2d_clawpatch_get_area(glob,patch);  
    clawpatch_vt->d2->fort_conservation_check(&mx, &my, &mbc, &mfields, &dx,&dy,
                                              area, rhs, error_data->rhs,
                                              error_data->c_kahan);
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);

    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);

    int intersects_bc[4];
    fclaw_physical_get_bc(glob,blockno,patchno,intersects_bc);

    double t = glob->curr_time;
    int cons_check = 1;

    POISSON_FORT_APPLY_BC(&blockno, &mx, &my, &mbc, &mfields, 
                         &xlower, &ylower, &dx,&dy,&t, intersects_bc,
                         mg_opt->boundary_conditions,rhs, mg_vt->fort_eval_bc,
                         &cons_check, error_data->boundary);

}


static
void poisson_time_header_ascii(fclaw_global_t* glob, int iframe)
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
void cb_poisson_output_ascii(fclaw_domain_t * domain,
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
    fclaw_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,&local_num, &level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    double *err;
    fclaw_clawpatch_elliptic_error_data(glob,patch,&err,&mfields);

    double *soln;
    fclaw_clawpatch_elliptic_soln_data(glob,patch,&soln,&mfields);

    char fname[BUFSIZ];
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);

    /* The fort routine is defined by a clawpack solver and handles 
       the layout of q in memory (i,j,m) or (m,i,j), etc */
    //fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    //FCLAW_ASSERT(clawpatch_vt->fort_output_ascii);

    POISSON_FORT_OUTPUT_ASCII(fname,&mx,&my,&mfields,&mbc,
                             &xlower,&ylower,&dx,&dy,rhs,
                             soln, err,
                             &global_num,&level,&blockno,
                             &glob->mpirank);
}


int poisson_tag4refinement(fclaw_global_t *glob,
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
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw_clawpatch_rhs_data(glob,this_patch,&rhs,&mfields);

    tag_patch = 0;
    fclaw_global_set_static(glob);
    clawpatch_vt->d2->fort_tag4refinement(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                          &blockno, rhs,&refine_threshold,
                                          &initflag,&tag_patch);
    fclaw_global_clear_static();
    return tag_patch;
}

static
int poisson_tag4coarsening(fclaw_global_t *glob,
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
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs[4];
    int mfields;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw_clawpatch_rhs_data(glob,&fine_patches[igrid],&rhs[igrid],&mfields);
    }

    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

    int tag_patch = 0;
    fclaw_global_set_static(glob);
    clawpatch_vt->d2->fort_tag4coarsening(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs[0],rhs[1],rhs[2],rhs[3],
                                      &coarsen_threshold,&initflag,&tag_patch);
    fclaw_global_clear_static();
    return tag_patch == 1;
}


void poisson_link_solvers(fclaw_global_t *glob)
{
    /* ForestClaw vtable */
    fclaw_vtable_t *fc_vt = fclaw_vt(glob);
    fc_vt->problem_setup = &poisson_problem_setup;  

    /* Patch : RHS function */
    fclaw_patch_vtable_t* patch_vt = fclaw_patch_vt(glob);
    patch_vt->rhs = poisson_rhs;          /* Overwrites default */
    patch_vt->initialize = poisson_rhs;   /* Get an initial refinement */

    /* Multigrid vtable */
    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);
    mg_vt->fort_rhs       = &POISSON_FORT_RHS;

    mg_vt->fort_beta      = &POISSON_FORT_BETA;
    
    mg_vt->fort_apply_bc = &POISSON_FORT_APPLY_BC;
    mg_vt->fort_eval_bc  = &POISSON_FORT_EVAL_BC;   // For non-homogeneous BCs

    /* Clawpatch : Compute the error */
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->compute_error = poisson_compute_error;
    clawpatch_vt->d2->fort_compute_patch_error = &POISSON_COMPUTE_ERROR;

    // tagging routines
    patch_vt->tag4refinement       = poisson_tag4refinement;
    patch_vt->tag4coarsening       = poisson_tag4coarsening;
#if 0
    clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;
#endif    

    // Setup routine
    patch_vt->setup = poisson_patch_setup;

    // Output routines

    fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    if (fclaw_opt->compute_error) 
    {
        clawpatch_vt->time_header_ascii = poisson_time_header_ascii;
        //clawpatch_vt->fort_header_ascii = &poisson_FORT_HEADER_ASCII;
        clawpatch_vt->cb_output_ascii = cb_poisson_output_ascii;        
    }

    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);

    if (mg_opt->patch_operator == USER_OPERATOR)
    {
        mg_vt->patch_operator = fc2d_thunderegg_fivepoint_solve;        
    }

    /* Diagnostics - most need to be modified because different number of fields are
       used for the elliptic problem then for the hyperbolic problem */
    clawpatch_vt->conservation_check = poisson_conservation_check;        

    fclaw_diagnostics_vtable_t *diag_vt = fclaw_diagnostics_vt(glob);
    diag_vt->patch_init_diagnostics     = poisson_diagnostics_initialize;
    diag_vt->patch_reset_diagnostics    = poisson_diagnostics_reset;
    diag_vt->patch_compute_diagnostics  = poisson_diagnostics_compute;
    diag_vt->patch_gather_diagnostics   = poisson_diagnostics_gather;
    diag_vt->patch_finalize_diagnostics = poisson_diagnostics_finalize;
}

