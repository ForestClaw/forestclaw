/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "swirlcons_user.h"

#include <fclaw_include_all.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>    /* Headers for user defined fortran files */

#include <fclaw_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>  /* headers for tag2refinement, tag4coarsening  */

#include "../all/advection_user_fort.h"  

static
void cb_swirl_output_ascii (fclaw_domain_t * domain,
                            fclaw_patch_t * this_patch,
                            int this_block_idx, int this_patch_idx,
                            void *user);


void swirlcons_link_solvers(fclaw_global_t *glob)
{
    fclaw2d_vtable_t                     *vt = fclaw2d_vt(glob);
    fclaw2d_patch_vtable_t         *patch_vt = fclaw2d_patch_vt(glob);
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt(glob);

    fc2d_clawpack46_options_t  *clawopt = fc2d_clawpack46_get_options(glob);
    const user_options_t*          user = swirlcons_get_options(glob);

    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    /* ForestClaw functions */
    vt->problem_setup = &swirlcons_problem_setup;  /* Version-independent */

    patch_vt->setup   = &swirlcons_patch_setup_manifold;

    clawpatch_vt->fort_compute_patch_error = &SWIRL46_COMPUTE_ERROR;

    clawopt->use_fwaves = 0;
    switch(user->mapping)
    {
        case 0:

            printf("Version without mapping does not work\n");
            exit(0);

            clawpack46_vt->fort_rpt2      = &RPT2CONS;
            clawpack46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE;
            switch(user->rp_solver)
            {
                case 1:  /* QS */
                    //clawpack46_vt->fort_rpn2      = &RPN2CONS_QS;
                case 2:  /* WD */
                    clawpack46_vt->fort_rpn2      = &RPN2CONS_WD; 
                case 3:  /* EC */
                    clawpack46_vt->fort_rpn2      = &RPN2CONS_EC;                 
                case 4:  /* FW */
                    clawopt->use_fwaves = 1;
                    clawpack46_vt->fort_rpn2      = RPN2CONS_FW;
            }
            // clawpack46_vt->fort_setaux = &CLAWPACK46_SETAUX;
            break;

        case 1: /* Cart */
        case 2: /* fivepatch */
        case 3: /* bilinear */
            clawpack46_vt->fort_rpt2      = &RPT2CONS_MANIFOLD;      
            clawpack46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;
            switch(user->rp_solver)
            {
                case 1:  /* QS */
                    //clawpack46_vt->fort_rpn2      = &RPN2CONS_QS_MANIFOLD; 
                case 2:  /* WD */
                    clawpack46_vt->fort_rpn2      = &RPN2CONS_WD; 
                case 3:  /* EC */
                    clawpack46_vt->fort_rpn2      = &RPN2CONS_EC_MANIFOLD;                 
                case 4:  /* FW */
                    clawopt->use_fwaves = 1;
                    clawpack46_vt->fort_rpn2      = RPN2CONS_FW_MANIFOLD; 
            }

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
                    clawpatch_vt->fort_header_ascii   = &SWIRL46_FORT_HEADER_ASCII;
                    clawpatch_vt->cb_output_ascii     = &cb_swirl_output_ascii;                
                }
            }

            break;
    }

    clawpack46_vt->fort_qinit     = CLAWPACK46_QINIT;
}

void swirlcons_problem_setup(fclaw_global_t* glob)
{
    const user_options_t* user = swirlcons_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->mapping,"\% mapping\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24d   %s",user->color_equation,"\% color_equation\n");
        fprintf(f,  "%-24d   %s",user->use_stream,"\% use_stream\n");
        fprintf(f,  "%-24.16f   %s",user->alpha,"\% alpha\n");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
 
    fclaw2d_domain_barrier (glob->domain);  /* redundant?  */
    SWIRLCONS_SETPROB();

#if 0
    SWIRL_SETPROB(&user->example,&user->mapping,&user->initial_condition,
                  &user->color_equation, &user->use_stream,
                  &user->alpha);
#endif
}


void swirlcons_patch_setup_manifold(fclaw_global_t *glob,
                                    fclaw_patch_t *patch,
                                    int blockno,
                                    int patchno)
{
    const user_options_t* user = swirlcons_get_options(glob);


    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_grid_data_2d(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);


    double *xp, *yp, *zp, *xd, *yd, *zd, *area;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *edgelengths, *curvature;
    fclaw_clawpatch_metric_scalar_2d(glob, patch, &area,&edgelengths,
                                    &curvature);

    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw_clawpatch_metric_vector_2d(glob,patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);

    switch(user->mapping)
    {
        case 0:  /* No map */
#if 0            
            int maxmx, maxmy;
            maxmx = mx;
            maxmy = my;
            CLAWPACK46_SETAUX(&maxmx, &maxmy, &mbc,&mx,&my,&xlower,&ylower,
                              &dx,&dy,&maux,aux);
#endif                              
            break;
        case 1: /* cart map */
        case 2: /* fivepatch */
        case 3: /* bilinear */
            SWIRL46_SETAUX(&mbc,&mx,&my,&xlower,&ylower,
                           &dx,&dy,&maux,aux,&blockno,
                           area, edgelengths,xnormals,ynormals,
                           surfnormals);
            break;
    }
}

static
void cb_swirl_output_ascii (fclaw_domain_t * domain,
                            fclaw_patch_t * patch,
                            int blockno, int patchno, 
                            void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_global_t *glob = (fclaw_global_t*) s->glob;

    int iframe = *((int *) s->user);
    double time = glob->curr_time;


    /* Get info not readily available to user */
    int local_patch_num, global_num, level;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num, 
                           &local_patch_num,&level);
    
    int mx,my,mbc;
    double xlower, ylower, dx, dy;
    fclaw_clawpatch_grid_data_2d(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);


    double* error = fclaw2d_clawpatch_get_error(glob,patch);
    double* soln = fclaw2d_clawpatch_get_exactsoln(glob,patch);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);


    /* Here, we pass in q and the error, so need special headers and files */
    SWIRL46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                            &xlower,&ylower,
                            &dx,&dy,
                            q,error,soln, &time, 
                            &local_patch_num,&level,
                            &blockno,
                            &glob->mpirank);
}









