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

#include "sphere_user.h"

#include <fclaw2d_include_all.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>    /* Headers for user defined fortran files */

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>  /* headers for tag2refinement, tag4coarsening  */

#include <fclaw2d_block.h>

static
void sphere_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = sphere_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->mapping,"\% mapping\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24.16f   %s",user->omega[0],"\% omega[0]\n");
        fprintf(f,  "%-24.16f   %s",user->omega[1],"\% omega[1]\n");
        fprintf(f,  "%-24.16f   %s",user->omega[2],"\% omega[2]\n");
        fprintf(f,  "%-24d   %s",user->refine_pattern,"\% refinement_pattern\n");
        fprintf(f,  "%-24.16f   %s",user->b_init,"\% b_init\n");
        fprintf(f,  "%-24.16f   %s",user->c_init,"\% c_init\n");
        fprintf(f,  "%-24.16f   %s",user->sharpness,"\% sharpness\n");
        fprintf(f,  "%-24d   %s",user->curvature_correction,"\% Curvature correction\n");
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    SPHERE_SETPROB();
}

void sphere_patch_setup_manifold(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    int blockno,
                                    int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xp, *yp, *zp, *xd, *yd, *zd, *area;
    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *edgelengths,*curvature;
    fclaw2d_clawpatch_metric_scalar(glob, this_patch,&area,&edgelengths,
                                    &curvature);

    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw2d_clawpatch_metric_vector(glob,this_patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    int maux;
    double *aux;
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    SPHERE_SETAUX(&blockno, &mx,&my,&mbc, &xlower,&ylower,
                  &dx,&dy, area, edgelengths,xp,yp,zp,
                  aux, &maux);
}


static
void sphere_b4step2(fclaw2d_global_t *glob,
                    fclaw2d_patch_t *this_patch,
                    int blockno,
                    int patchno,
                    double t, double dt)

{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw2d_clawpatch_metric_vector(glob,this_patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    SPHERE_SET_VELOCITIES(&blockno, &mx, &my, &mbc,
                          &dx, &dy, &xlower, &ylower,
                          &t, xnormals,ynormals, surfnormals,
                          aux,&maux);

}

static
int sphere_tag4refinement(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *patch,
                          int blockno, int patchno,
                          int initflag)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double refine_threshold = fclaw_opt->refine_threshold;

    int tag_patch = 1;
    if (refine_threshold < 0) 
    {
        /* Always refine if refine-threshold < 0 */
        tag_patch = 1;
    }
    else
    {
        /* Pass in time to tagging routine */
        tag_patch = 0;  
        double t = glob->curr_time;
        SPHERE_TAG4REFINEMENT(&mx,&my,&mbc,&meqn, &xlower,&ylower,&dx,&dy,
                              &t, &blockno, q,&refine_threshold,
                              &initflag,&tag_patch);
    }
    return tag_patch;
}


static
int sphere_tag4coarsening(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *fine_patches,
                          int blockno,
                          int patchno,
                          int initflag)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_patch_t *patch0 = &fine_patches[0];
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q[4];
    int meqn;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
    }

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int tag_patch = 0;
    if (coarsen_threshold >= 0) 
    {
        tag_patch = 0;
        double t = glob->curr_time;
        SPHERE_TAG4COARSENING(&mx,&my,&mbc,&meqn,
                              &xlower,&ylower,&dx,&dy,&t,
                              &blockno, q[0],q[1],q[2],q[3],
                              &coarsen_threshold,&initflag,&tag_patch);
    }
    return tag_patch == 1;
}



static
void cb_sphere_output_ascii (fclaw2d_domain_t * domain,
                            fclaw2d_patch_t * this_patch,
                            int this_block_idx, int this_patch_idx,
                            void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t      *glob = (fclaw2d_global_t*) s->glob;

    int iframe = *((int *) s->user);
    double time = glob->curr_time;


    /* Get info not readily available to user */
    int local_num, global_num, level;
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &global_num, &local_num,&level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *q;
    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    double* error = fclaw2d_clawpatch_get_error(glob,this_patch);
    double* soln = fclaw2d_clawpatch_get_exactsoln(glob,this_patch);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);


    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);


    /* Here, we pass in q and the error, so need special headers and files */
    SPHERE_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,&maux,
                            &xlower,&ylower,
                            &dx,&dy,
                            q,error,soln, aux, 
                            &time, &global_num,&level,
                            &this_block_idx,
                            &fclaw_opt->compute_error, 
                            &glob->mpirank);
}



void sphere_link_solvers(fclaw2d_global_t *glob)
{
    /* ForestClaw core functions */
    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &sphere_problem_setup;  /* Version-independent */

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
    patch_vt->setup   = &sphere_patch_setup_manifold;
    patch_vt->tag4refinement = sphere_tag4refinement;
    patch_vt->tag4coarsening = sphere_tag4coarsening;

    fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt(glob);
    clawpack46_vt->b4step2        = sphere_b4step2;
    clawpack46_vt->fort_qinit     = CLAWPACK46_QINIT;
    clawpack46_vt->fort_rpn2      = RPN2CONS_FW_MANIFOLD; 
    clawpack46_vt->fort_rpt2      = &RPT2CONS_MANIFOLD;      
    clawpack46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;

    /* Include error in output files */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    clawpatch_vt->fort_header_ascii   = &SPHERE_FORT_HEADER_ASCII;
    clawpatch_vt->cb_output_ascii     = &cb_sphere_output_ascii;                

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        clawpatch_vt->fort_compute_patch_error = &SPHERE_COMPUTE_ERROR;
    }
 }










