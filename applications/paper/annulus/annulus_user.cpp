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

#include "annulus_user.h"

#include <fclaw2d_include_all.h>

#include "fclaw2d_clawpatch.h"

#include <fc2d_clawpack46.h>

#include "fclaw2d_clawpatch_options.h"
#include <fc2d_clawpack46_options.h>


#include "clawpack46_advection_user_fort.h"


static
void annulus_problem_setup(fclaw2d_global_t *glob)
{
    const user_options_t *user = annulus_get_options(glob);

    if (glob->mpirank == 0)
    {

        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,       "\% example\n");    
        fprintf(f,  "%-24d   %s",user->mapping,       "\% mapping\n");   
        fprintf(f,  "%-24d   %s",user->initchoice,    "\% initchoice\n");    
        fprintf(f,"%-24.4f   %s",user->revs_per_s,    "\% revs_per_s\n");    
        fprintf(f,"%-24.4f   %s",user->twist,         "\% twist\n");    
        fprintf(f,"%-24.4f   %s",user->cart_speed,    "\% cart_speed\n");    
        fprintf(f,  "%-24d   %s",user->color_equation,"\% color_equation\n");    
        fprintf(f,  "%-24d   %s",user->use_stream,    "\% use_stream\n");    
        fprintf(f,"%-24.4f   %s",user->beta,          "\% beta\n");    
        fprintf(f,"%-12.4f%-12.4f    %s",user->theta[0],user->theta[1],"\% beta\n");    
        fprintf(f,  "%-24d   %s",user->refine_pattern,"\% refine_pattern\n");    
        fprintf(f,"%-24.4f   %s",user->init_radius,   "\% init_radius\n");    
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    SETPROB_ANNULUS();  /* Reads file created above */
 
}


static
void annulus_patch_setup(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int blockno,
                         int patchno)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *curvature, *surfnormals;
    double *edgelengths;
    double *xnormals, *ynormals, *xtangents, *ytangents;
    const user_options_t* user = annulus_get_options(glob);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    fclaw2d_clawpatch_metric_scalar(glob, this_patch,&area,&edgelengths,
                                    &curvature);

    fclaw2d_clawpatch_metric_vector(glob,this_patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    if (user->claw_version == 4)
    {
        ANNULUS46_SETAUX(&mbc,&mx,&my,&xlower,&ylower,
                         &dx,&dy,&maux,aux,&blockno,
                         area,xd,yd,zd,
                         edgelengths,xnormals,ynormals,
                         xtangents, ytangents, surfnormals);
    }
#if 0
    else if(user->claw_version == 5)
    {
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                              &maux,aux,&blockno,xd,yd,zd,area);
    }
#endif    
}


static
void cb_annulus_output_ascii (fclaw2d_domain_t * domain,
                              fclaw2d_patch_t * this_patch,
                              int this_block_idx, int this_patch_idx,
                              void *user)
{
    int patch_num;
    int level;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy, time;
    double *q, *error, *soln;
    int iframe;

    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t      *glob = (fclaw2d_global_t*) s->glob;

    //fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    const fclaw_options_t         *fclaw_opt = fclaw2d_get_options(glob);


    iframe = *((int *) s->user);

    time = glob->curr_time;


    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &patch_num,&level);
    
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    error = fclaw2d_clawpatch_get_error(glob,this_patch);
    soln = fclaw2d_clawpatch_get_exactsoln(glob,this_patch);

    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);


    /* Here, we pass in q and the error, so need special headers and files */
    ANNULUS46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                              &xlower,&ylower,
                              &dx,&dy,
                              q,error,soln, &time, 
                              &patch_num,&level,
                              &this_block_idx,
                              &glob->mpirank);
}


void annulus_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t           *vt            = fclaw2d_vt();
    fclaw2d_patch_vtable_t     *patch_vt      = fclaw2d_patch_vt();
    fclaw2d_clawpatch_vtable_t *clawpatch_vt  = fclaw2d_clawpatch_vt();
    //fc2d_clawpack46_vtable_t   *clawpack46_vt = fc2d_clawpack46_vt();

    fclaw_options_t          *fclaw_opt = fclaw2d_get_options(glob);
    const user_options_t          *user = annulus_get_options(glob);

    vt->problem_setup  = &annulus_problem_setup;
    patch_vt->setup    = &annulus_patch_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_options_t  *clawopt     = fc2d_clawpack46_get_options(glob);
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();

        clawpack46_vt->fort_qinit   = CLAWPACK46_QINIT;
#if 0 
        /* Doesn't really work with transverse solvers */
        clawpack46_vt->fort_bc2     = CLAWPACK46_BC2;  /* Replace default version */
#endif        

        clawpatch_vt->fort_compute_patch_error = &ANNULUS46_COMPUTE_ERROR;
        clawpatch_vt->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

        clawopt->use_fwaves = 1;
        clawpack46_vt->fort_rpn2      = &RPN2CONS_FW_MANIFOLD;   
        clawpack46_vt->fort_rpt2      = &RPT2CONS_MANIFOLD;      
        //clawpack46_vt->fort_rpt2      = &ANNULUS46_RPT2ADV_MANIFOLD;      
        clawpack46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;


        if (fclaw_opt->compute_error)
        {
            clawpatch_vt->fort_header_ascii   = &ANNULUS46_FORT_HEADER_ASCII;
            clawpatch_vt->cb_output_ascii     = &cb_annulus_output_ascii;                
        }
    }
#if 0    
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt();
        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        claw5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        claw5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;
    }
#endif    
}


