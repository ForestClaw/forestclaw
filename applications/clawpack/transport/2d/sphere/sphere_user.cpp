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
        fprintf(f,  "%-24.16f   %s",user->kappa,"\% kappa\n");
        fprintf(f,  "%-24.16f   %s",user->period,"\% period\n");
        fprintf(f,  "%-24.16f   %s",user->omega[0],"\% omega[0]\n");
        fprintf(f,  "%-24.16f   %s",user->omega[1],"\% omega[1]\n");
        fprintf(f,  "%-24.16f   %s",user->omega[2],"\% omega[2]\n");
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
    double xlower, ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double refine_threshold = fclaw_opt->refine_threshold;

    double *area, *edgelengths, *curvature;
    fclaw2d_clawpatch_metric_scalar(glob, patch,&area,&edgelengths,
                                    &curvature);

    int tag_patch = 0;
#if 0
    //int *corner_count = fclaw2d_patch_block_corner_count(glob,patch);

    int intersects_block[4];
    fclaw2d_block_get_block_boundary(glob, patch, intersects_block);
    int block_edge = 0;
    for(int k = 0; k < 4; k++)
    {
        if (intersects_block[k] == 1)
        {
            block_edge = 1;
            break;
        }
    }
    if (block_edge)
    {
        tag_patch = 1;
    }
#endif
    SPHERE_TAG4REFINEMENT(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                          &blockno, q, curvature, &refine_threshold,
                          &initflag,&tag_patch);
    return tag_patch;
}

static
int sphere_tag4coarsening(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *fine_patches,
                          int blockno,
                          int patchno)
{
    //fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();


    int tag_patch,igrid;
    fclaw2d_patch_t *patch0;

    patch0 = &fine_patches[0];

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *c[4], *q[4];
    for (igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);

        double *area, *edgelengths;
        fclaw2d_clawpatch_metric_scalar(glob, &fine_patches[igrid],
                                        &area,&edgelengths,
                                        &c[igrid]);

    }

    tag_patch = 0;

#if 0    
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;
    SPHERE_TAG4COARSENING(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                          &blockno, 
                          q[0],q[1],q[2],q[3],
                          c[0],c[1],c[2],c[3],
                          &coarsen_threshold,&tag_patch);
#endif                          
    return tag_patch == 1;
}



static
void cb_sphere_output_ascii (fclaw2d_domain_t * domain,
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
    SPHERE_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                            &xlower,&ylower,
                            &dx,&dy,
                            q,error,soln, &time, 
                            &patch_num,&level,
                            &this_block_idx,
                            &glob->mpirank);
}



void sphere_link_solvers(fclaw2d_global_t *glob)
{
    /* ForestClaw core functions */
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    vt->problem_setup = &sphere_problem_setup;  /* Version-independent */

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    patch_vt->setup   = &sphere_patch_setup_manifold;
    patch_vt->tag4refinement = sphere_tag4refinement;
    patch_vt->tag4coarsening = sphere_tag4coarsening;

    fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt();
    clawpack46_vt->b4step2        = sphere_b4step2;
    clawpack46_vt->fort_qinit     = CLAWPACK46_QINIT;
    clawpack46_vt->fort_rpn2      = RPN2CONS_FW_MANIFOLD; 
    clawpack46_vt->fort_rpt2      = &RPT2CONS_MANIFOLD;      
    clawpack46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;

    /* Clawpatch functions */
    //const user_options_t* user = sphere_get_options(glob);
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    // clawpatch_vt->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
    //clawpatch_vt->fort_tag4coarsening = &SPHERE_TAG4COARSENING;        

    /* Include error in output files */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        clawpatch_vt->fort_compute_patch_error = &SPHERE_COMPUTE_ERROR;
        clawpatch_vt->fort_header_ascii   = &SPHERE_FORT_HEADER_ASCII;
        clawpatch_vt->cb_output_ascii     = &cb_sphere_output_ascii;                
    }
 }










