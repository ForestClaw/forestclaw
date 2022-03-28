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

#include "square_user.h"

#include <fclaw2d_include_all.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>    /* Headers for user defined fortran files */

#include <fc2d_clawpack5.h>
#include <fc2d_clawpack5_options.h>
#include <clawpack5_user_fort.h>    /* Headers for user defined fortran files */

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>  /* headers for tag2refinement, tag4coarsening  */

static
void square_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = square_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->mapping,"\% mapping\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,"\% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->center[0],"\% x0\n");
        fprintf(f,  "%-24.6f   %s",user->center[1],"\% y0\n");
        fprintf(f,  "%-24.6f   %s",user->velocity[0],"\% u\n");
        fprintf(f,  "%-24.6f   %s",user->velocity[1],"\% v\n");
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    SQUARE_SETPROB();


}

void square_patch_setup_manifold(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *patch,
                                    int blockno,
                                    int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xp, *yp, *zp, *xd, *yd, *zd, *area;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *edgelengths,*curvature;
    fclaw2d_clawpatch_metric_scalar(glob, patch,&area,&edgelengths,
                                    &curvature);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);    

    const user_options_t* user = square_get_options(glob);
    if (user->claw_version == 4) 
    {
        SQUARE46_SETAUX(&blockno, &mx,&my,&mbc, &xlower,&ylower,
            &dx,&dy, area, edgelengths,xp,yp,zp,
            aux, &maux);
    }
    else if (user->claw_version == 5)
    {
        SQUARE5_SETAUX(&blockno, &mx,&my,&mbc, &xlower,&ylower,
            &dx,&dy, area, edgelengths,xp,yp,zp,
            aux, &maux);
    }


    /* Assume that velocities don't depend on t */
    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw2d_clawpatch_metric_vector(glob,patch, &xnormals, &ynormals, 
                                    &xtangents, &ytangents, &surfnormals);

    double t = 0; /* Not used */
    if (user->claw_version == 4)
    {
        SQUARE46_SET_VELOCITIES(&blockno, &mx, &my, &mbc,
                                &dx, &dy, &xlower, &ylower,
                                &t, xnormals,ynormals, surfnormals,
                                aux,&maux);
    }
    else
    {
        SQUARE5_SET_VELOCITIES(&blockno, &mx, &my, &mbc,
                                &dx, &dy, &xlower, &ylower,
                                &t, xnormals,ynormals, surfnormals,
                                aux,&maux);        
    }
}

static
void cb_square_output_ascii (fclaw2d_domain_t * domain,
                            fclaw2d_patch_t * patch,
                            int blockno, int patchno,
                            void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t  *glob = (fclaw2d_global_t*) s->glob;
    const fclaw_options_t  *fclaw_opt = fclaw2d_get_options(glob);

    int iframe = *((int *) s->user);
    double time = glob->curr_time;

    /* Get info not readily available to user */
    int level, patch_num, global_num;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num, &patch_num,&level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);
    double* error = fclaw2d_clawpatch_get_error(glob,patch);
    double* soln = fclaw2d_clawpatch_get_exactsoln(glob,patch);

    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);


    /* Here, we pass in q and the error, so need special headers and files */
    const user_options_t* user_opt = square_get_options(glob);
    if (user_opt->claw_version == 4)
    {
        SQUARE46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
            &xlower,&ylower, &dx,&dy, q,error,soln, &time,
            &patch_num,&level, &blockno, &glob->mpirank);
    }
    else if (user_opt->claw_version == 5)
    {
        SQUARE5_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
            &xlower,&ylower, &dx,&dy, q,error,soln, &time,
            &patch_num,&level, &blockno, &glob->mpirank);        
    }
}



void square_link_solvers(fclaw2d_global_t *glob)
{
    /* ForestClaw core functions */
    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &square_problem_setup;  /* Version-independent */

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
    patch_vt->setup = &square_patch_setup_manifold;

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);

    const user_options_t* user = square_get_options(glob);
    if (user->claw_version == 4) 
    {
        fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt();
        clawpack46_vt->fort_qinit  = &CLAWPACK46_QINIT;
        
        /* Be careful : Signatures for rpn2fw, rpt2fw not the same as rpn2 and rpt2fw. */
        clawpack46_vt->fort_rpn2fw   = &CLAWPACK46_RPN2FW_MANIFOLD; 
        clawpack46_vt->fort_rpt2fw   = &CLAWPACK46_RPT2FW_MANIFOLD;      
        clawpack46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;

        /* Clawpatch functions */
        //const user_options_t* user = square_get_options(glob);
        clawpatch_vt->fort_tag4refinement = &SQUARE46_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &SQUARE46_TAG4COARSENING;       
    } 
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t  *clawpack5_vt = fc2d_clawpack5_vt();
        clawpack5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        clawpack5_vt->fort_rpn2      = &CLAWPACK5_RPN2FW_MANIFOLD; 
        clawpack5_vt->fort_rpt2      = &CLAWPACK5_RPT2_MANIFOLD;      
        clawpack5_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;

        /* Clawpatch functions */
        clawpatch_vt->fort_tag4refinement = &SQUARE5_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &SQUARE5_TAG4COARSENING;     
    }

#if 0
    if (0)
    {        
        clawpatch_vt->fort_interpolate2fine      = SQUARE_FORT_INTERPOLATE2FINE;
        clawpatch_vt->fort_interpolate_face      = SQUARE_FORT_INTERPOLATE_FACE;
        clawpatch_vt->fort_interpolate_corner    = SQUARE_FORT_INTERPOLATE_CORNER;
    }
#endif

    /* Include error in output files */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        if (user->claw_version == 4)
        {
            clawpatch_vt->fort_compute_patch_error = &SQUARE46_COMPUTE_ERROR;
            clawpatch_vt->fort_header_ascii   = &SQUARE46_FORT_HEADER_ASCII;            
        }
        else if (user->claw_version == 5)
        {
            clawpatch_vt->fort_compute_patch_error = &SQUARE5_COMPUTE_ERROR;
            clawpatch_vt->fort_header_ascii   = &SQUARE5_FORT_HEADER_ASCII;                        
        }
        clawpatch_vt->cb_output_ascii     = &cb_square_output_ascii;                
    }

 }










