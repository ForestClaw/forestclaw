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

#include "cylinder_user.h"

#include <fclaw_include_all.h>

#include <fclaw_clawpatch.h>
#include <fclaw2d_metric.h>
#include <fclaw_clawpatch_options.h>

/* Two versions of Clawpack */
#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack46.h>

#include "clawpack46_advection_user_fort.h"


#if 0
void cylinder_compute_area(fclaw_global_t *glob,
                           fclaw_patch_t *this_patch,
                           int blockno,
                           int patchno);
#endif                           


static
void cylinder_problem_setup(fclaw_global_t *glob)
{
    const user_options_t* user = cylinder_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24d   %s",user->refine_pattern,"\% refine_pattern\n");
        fprintf(f,  "%-24.16f   %s",user->R,"\% R\n");
        fprintf(f,  "%-24.16f   %s",user->H,"\% H\n");
        fprintf(f,  "%-24.16f   %s",user->xc0,"\% xc0\n");
        fprintf(f,  "%-24.16f   %s",user->yc0,"\% yc0\n");
        fprintf(f,  "%-24.16f   %s",user->r0,"\% r0\n");
        fprintf(f,  "%-24.16f   %s",user->revs_per_s,"\% revs_per_second\n");
        fprintf(f,  "%-24.16f   %s",user->v_speed,"\% v_speed\n");
        fprintf(f,  "%-24d      %s",user->exact_metric,"\% exact_metric\n");
        fprintf(f,  "%-24d      %s",user->mapping,"\% mapping\n");
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    CYLINDER_SETPROB();        
}


static
void cylinder_patch_setup(fclaw_global_t *glob,
                       fclaw_patch_t *patch,
                       int blockno, int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_grid_data_2d(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xp, *yp, *zp, *xd, *yd, *zd, *area;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *edgelengths,*curvature;
    fclaw2d_clawpatch_metric_scalar(glob, patch,&area,&edgelengths,
                                    &curvature);

    int maux;
    double *aux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);
    CYLINDER_SETAUX(&blockno, &mx,&my,&mbc, &xlower,&ylower,
                  &dx,&dy, area, edgelengths,xp,yp,zp,
                  aux, &maux);


    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw2d_clawpatch_metric_vector(glob,patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    double t = 0; /* Not used */
    CYLINDER_SET_VELOCITIES(&blockno, &mx, &my, &mbc,
                         &dx, &dy, &xlower, &ylower,
                         &t, xnormals,ynormals, surfnormals,
                         aux,&maux);
}


static
void cb_cylinder_output_ascii (fclaw_domain_t * domain,
                            fclaw_patch_t * this_patch,
                            int blockno, int patchno,
                            void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_global_t      *glob = (fclaw_global_t*) s->glob;

    const user_options_t *user_opt =  cylinder_get_options(glob);

    int iframe = *((int *) s->user);
    double time = glob->curr_time;


    /* Get info not readily available to user */
    int level,patch_num, global_num;
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           blockno,patchno,
                           &global_num, &patch_num,&level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_grid_data_2d(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double* q;
    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    double* error = fclaw2d_clawpatch_get_error(glob,this_patch);
    double* soln = fclaw2d_clawpatch_get_exactsoln(glob,this_patch);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);


    /* Here, we pass in q and the error, so need special headers and files */
    if (user_opt->claw_version == 4)
    {
        CYLINDER46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                                &xlower,&ylower,
                                &dx,&dy,
                                q,error,soln, &time, 
                                &patch_num,&level,
                                &blockno,
                                &glob->mpirank);
    }
#if 0    
    else if (user_opt->claw_version == 5)
    {
        CYLINDER5_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,&xlower,&ylower,
                               &dx,&dy,q,error,
                               &patch_num,&level,&blockno,
                               &glob->mpirank);
    }
#endif    
}

static
void cylinder_compute_area(fclaw_global_t *glob,
                           fclaw_patch_t *patch,
                           int blockno,
                           int patchno)
{
    const fclaw_options_t* gparms = fclaw_get_options(glob);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_metric_patch_grid_data(glob,patch,&mx,&my,&mbc,
                                   &xlower,&ylower,&dx,&dy);

    double *area = fclaw2d_metric_patch_get_area(patch);

    int maxlevel = gparms->maxlevel;
    int level = patch->level;

    const user_options_t* user = cylinder_get_options(glob);
    if (user->mapping == 0)
        CYLINDER_COMPUTE_AREA(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                              &blockno, &maxlevel, &level, area);
    else
        LATLONG_COMPUTE_AREA(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                              &blockno, &maxlevel, &level, area);
}


static
void cylinder_compute_tensors(fclaw_global_t *glob,
                              fclaw_patch_t *patch,
                              int blockno,
                              int patchno)
{
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt(glob);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *xnormals, *ynormals;
    double *xtangents, *ytangents;
    double *edgelengths;
    double *surfnormals, *curvature;
    double *area;

    fclaw2d_metric_patch_grid_data(glob,patch,&mx,&my,&mbc,
                                   &xlower,&ylower,&dx,&dy);

    fclaw2d_metric_patch_mesh_data(glob,patch,
                                &xp,&yp,&zp,&xd,&yd,&zd,&area);

    fclaw2d_metric_patch_mesh_data2(glob,patch,
                                 &xnormals,&ynormals,
                                 &xtangents,&ytangents,
                                 &surfnormals,&edgelengths,
                                 &curvature);

    /* The user could set these to NULL to avoid doing these computations ... */

    const user_options_t* user = cylinder_get_options(glob);
    if (user->mapping == 0) 
    {
        CYLINDER_COMPUTE_NORMALS(&mx,&my,&mbc,xp,yp,zp,xd,yd,zd,
                                     xnormals,ynormals);

        int level=patch->level;
        CYLINDER_COMPUTE_TANGENTS(&mx,&my,&mbc,&dx, &dy, &level,
                                  xd,yd,zd,xtangents,ytangents,
                                  edgelengths);
    }       
    else
    {
        LATLONG_COMPUTE_NORMALS(&mx,&my,&mbc,xp,yp,zp,xd,yd,zd,
                                     xnormals,ynormals);

        int level=patch->level;
        LATLONG_COMPUTE_TANGENTS(&mx,&my,&mbc,&dx, &dy, &level,
                                  xd,yd,zd,xtangents,ytangents,
                                  edgelengths);

    }
    if (metric_vt->fort_compute_surf_normals != NULL)
    {
        metric_vt->fort_compute_surf_normals(&mx,&my,&mbc,xnormals,ynormals,edgelengths,
                                             curvature, surfnormals, area);
    }
}


void cylinder_link_solvers(fclaw_global_t *glob)
{

    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &cylinder_problem_setup;  /* Version-independent */

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
    patch_vt->setup   = &cylinder_patch_setup;

    fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);
    claw46_vt->fort_qinit = &CLAWPACK46_QINIT;
    claw46_vt->fort_rpn2  = RPN2CONS_FW_MANIFOLD; 
    claw46_vt->fort_rpt2  = &RPT2CONS_MANIFOLD;  
    claw46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->fort_tag4refinement = &CYLINDER_TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &CYLINDER_TAG4COARSENING;


    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt(glob);
    /* Area and edge lengths can be computed analytically  */
    metric_vt->compute_area          = cylinder_compute_area;
    // metric_vt->compute_area_ghost    = fclaw2d_metric_compute_area_ghost_default;
    metric_vt->compute_tensors       = cylinder_compute_tensors;

    /* Fortran files */
    // metric_vt->fort_compute_normals       = &FCLAW2D_FORT_COMPUTE_NORMALS;
    // metric_vt->fort_compute_tangents      = &FCLAW2D_FORT_COMPUTE_TANGENTS;
    // metric_vt->fort_compute_surf_normals  = &FCLAW2D_FORT_COMPUTE_SURF_NORMALS;



    /* Include error in output files */
    const fclaw_options_t  *fclaw_opt = fclaw_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        clawpatch_vt->fort_compute_patch_error = &CYLINDER46_COMPUTE_ERROR;
        clawpatch_vt->fort_header_ascii   = &CYLINDER46_FORT_HEADER_ASCII;
        clawpatch_vt->cb_output_ascii     = &cb_cylinder_output_ascii;                
    }

    /* Solve conservative equation using cell-centered velocity */

}







