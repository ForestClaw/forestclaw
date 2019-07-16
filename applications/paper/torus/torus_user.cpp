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

#include "torus_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

/* Two versions of Clawpack */
#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack46.h>

#include "clawpack46_advection_user_fort.h"


static
void torus_problem_setup(fclaw2d_global_t *glob)
{
    const user_options_t* user = torus_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->mapping,"\% mapping\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,"\% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->beta,"\% beta\n");
        fprintf(f,  "%-24.6f   %s",user->revs_per_s,"\% revs_per_second\n");
        fprintf(f,  "%-24d   %s",user->color_equation,"\% color_equation\n");    
        fprintf(f,  "%-24d   %s",user->use_stream,"\% use_stream\n");    
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    TORUS_SETPROB();
}


static
void torus_patch_setup(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    const user_options_t* user = torus_get_options(glob);

    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*edgelengths,*area, *curvature;
    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

#if 0
    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);
#endif                                  

    fclaw2d_clawpatch_metric_scalar(glob, this_patch,&area,&edgelengths,
                                    &curvature);

    fclaw2d_clawpatch_metric_vector(glob,this_patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    if (user->claw_version == 4)
    {
        /* Handles both non-conservative (ex 1-2) and conservative (ex 3-4) forms */
        FCLAW_ASSERT(maux == 7);
        TORUS46_SETAUX(&mbc,&mx,&my,&xlower,&ylower,
                       &dx,&dy,&maux,aux,&this_block_idx,
                       area, edgelengths,xnormals,ynormals, 
                       surfnormals);
    }
#if 0    
    else if (user->claw_version == 5)
    {
        /* Only works for advection in non-conservative form */
        fc2d_clawpack5_setaux(glob,this_patch,this_block_idx,this_patch_idx);
        fc2d_clawpack5_set_capacity(glob,this_patch,this_block_idx,this_patch_idx);
    }
#endif    
}

#if 0
void torus_patch_setup(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    const user_options_t* user = torus_get_options(glob);

    /* We use the non-mapped setaux routine, but need to set capacity */
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_setaux(glob,this_patch,this_block_idx,this_patch_idx);
        fc2d_clawpack46_set_capacity(glob,this_patch,this_block_idx,this_patch_idx);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_setaux(glob,this_patch,this_block_idx,this_patch_idx);
        fc2d_clawpack5_set_capacity(glob,this_patch,this_block_idx,this_patch_idx);
    }
}
#endif

static
void cb_torus_output_ascii (fclaw2d_domain_t * domain,
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
    const user_options_t           *user_opt =  torus_get_options(glob);


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
    if (user_opt->claw_version == 4)
    {
        TORUS46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                                &xlower,&ylower,
                                &dx,&dy,
                                q,error,soln, &time, 
                                &patch_num,&level,
                                &this_block_idx,
                                &glob->mpirank);
    }
#if 0    
    else if (user_opt->claw_version == 5)
    {
        TORUS5_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,&xlower,&ylower,
                               &dx,&dy,q,error,
                               &patch_num,&level,&this_block_idx,
                               &glob->mpirank);
    }
#endif    
}


void torus_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

    const fclaw_options_t         *fclaw_opt = fclaw2d_get_options(glob);

    vt->problem_setup = &torus_problem_setup;  /* Version-independent */
    patch_vt->setup   = &torus_patch_setup;

    const user_options_t *user =  torus_get_options(glob);
    clawpatch_vt->fort_compute_patch_error = &TORUS46_COMPUTE_ERROR;

    fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();
    fc2d_clawpack46_options_t *claw46_opt = fc2d_clawpack46_get_options(glob);

    claw46_vt->fort_qinit = &CLAWPACK46_QINIT;

    if (user->initial_condition == 0)
    {
            /* Initial condition is 0/1 field; Use absolute jump in solution */            
    }
    else
    {
        /* Smooth initial condition for accuracy problem : 
           We should use a divided differences for tagging */
        clawpatch_vt->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

        /* Include error in output files */
        if (fclaw_opt->compute_error)
        {
            clawpatch_vt->fort_header_ascii   = &TORUS46_FORT_HEADER_ASCII;
            clawpatch_vt->cb_output_ascii     = &cb_torus_output_ascii;                
        }

        if(user->color_equation)
        {
            claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV_MANIFOLD;
            claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV_MANIFOLD;

            /* Flux function (set to zero for non-conservative fix) */
            claw46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_ZERO;
        }
        else
        {
            /* Solve conservative equation using cell-centered velocity */
            /* Fwaves give best accuracy */
            claw46_opt->use_fwaves = 1;
            claw46_vt->fort_rpn2      = RPN2CONS_FW_MANIFOLD; 

            /* Transverse solver : Conservative form */
            claw46_vt->fort_rpt2      = &RPT2CONS_MANIFOLD;  

            /* Flux function (for conservative fix) */
            claw46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;
        }
    }
}



