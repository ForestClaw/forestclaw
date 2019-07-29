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

#include "swirlcons_user.h"

#include <fclaw2d_include_all.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>    /* Headers for user defined fortran files */

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>  /* headers for tag2refinement, tag4coarsening  */

#include "../all/advection_user_fort.h"  

static
void cb_swirl_output_ascii (fclaw2d_domain_t * domain,
                            fclaw2d_patch_t * this_patch,
                            int this_block_idx, int this_patch_idx,
                            void *user);


void swirlcons_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t                     *vt = fclaw2d_vt();
    fclaw2d_patch_vtable_t         *patch_vt = fclaw2d_patch_vt();
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt();

    fc2d_clawpack46_options_t  *clawopt = fc2d_clawpack46_get_options(glob);
    const user_options_t*          user = swirlcons_get_options(glob);

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

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

void swirlcons_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = swirlcons_get_options(glob);

    SWIRL_SETPROB(&user->example,&user->mapping,&user->initial_condition,
                  &user->color_equation, &user->use_stream,
                  &user->alpha);
}


void swirlcons_patch_setup_manifold(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    int blockno,
                                    int patchno)
{
    const user_options_t* user = swirlcons_get_options(glob);

    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*edgelengths,*area, *curvature;
    double *xp, *yp, *zp, *xd, *yd, *zd;
    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_scalar(glob, this_patch,&area,&edgelengths,
                                    &curvature);

    fclaw2d_clawpatch_metric_vector(glob,this_patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

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
                           area,
                           edgelengths,xnormals,ynormals,surfnormals);
            break;
    }
}

static
void cb_swirl_output_ascii (fclaw2d_domain_t * domain,
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
    SWIRL46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                            &xlower,&ylower,
                            &dx,&dy,
                            q,error,soln, &time, 
                            &patch_num,&level,
                            &this_block_idx,
                            &glob->mpirank);
}









