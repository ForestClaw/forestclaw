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

#include <fclaw2d_global.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw_math.h>
#include <fclaw2d_metric_default_fort.h>

void fclaw2d_metric_compute_area(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t* this_patch,
                                 int blockno,
                                 int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    int level, maxlevel,refratio;

    const amr_options_t* gparms = glob->gparms;
    level = this_patch->level;
    maxlevel = gparms->maxlevel;
    refratio = gparms->refratio;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *area = fclaw2d_clawpatch_get_area(glob,this_patch);

    /* Could make this a virtual function, but what is the signature?
       vt.fort_compute_area(...) */

    int m = pow_int(refratio,maxlevel-level);
    double *quadstore = FCLAW_ALLOC(double,3*(m+1)*(m+1));

    int ghost_only = 0;
    FCLAW2D_FORT_COMPUTE_AREA(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                              &blockno, area, &m, quadstore, &ghost_only);

    FCLAW_FREE(quadstore);
}

void fclaw2d_metric_compute_area_exact(fclaw2d_global_t *glob,
                                       fclaw2d_patch_t *this_patch,
                                       int blockno,
                                       int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *area = fclaw2d_clawpatch_get_area(glob,this_patch);

    int size = (2*(mbc+1) + mx)*(2*(mbc+1) + my);
    double *favg = FCLAW_ALLOC(double,size);

    int ghost_only = 0;
    int compute_avg = 0;
    fclaw2d_fort_aux_func_t f = &ONE;
    FCLAW2D_FORT_INTEGRATE_EXACT(&mx,&my,&mbc,&dx,&dy,
                                 &xlower, &ylower, &blockno,
                                 area, &f, favg, &compute_avg,
                                 &ghost_only);
    FCLAW_FREE(favg);
}


void fclaw2d_metric_area_set_ghost(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch,
                                   int blockno,
                                   int patchno)
{
    int mx,my, mbc;
    double xlower,ylower,dx,dy;
    double *area;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(glob,this_patch);

    /* Set area in ghost cells not set above */
    const amr_options_t* gparms = glob->gparms;
    int level = this_patch->level;
    int maxlevel = gparms->maxlevel;
    int refratio = gparms->refratio;
    int ghost_only = 1;

    int m = pow_int(refratio,maxlevel-level);
    double *quadstore = FCLAW_ALLOC(double,3*(m+1)*(m+1));

    FCLAW2D_FORT_COMPUTE_AREA(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                              &blockno, area, &m, quadstore,
                              &ghost_only);
    FCLAW_FREE(quadstore);
}

void fclaw2d_metric_area_set_ghost_exact(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch,
                                         int blockno,
                                         int patchno)
{
    int mx,my, mbc;
    double xlower,ylower,dx,dy;
    double *area;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(glob,this_patch);

    /* Not needed, but seems to be bad form to send in a NULL */
    int size = (2*(mbc+1) + mx)*(2*(mbc+1) + my);
    double *favg = FCLAW_ALLOC(double,size);

    /* Set area in ghost cells not set above */
    int ghost_only = 1;
    int compute_avg = 0;
    fclaw2d_fort_aux_func_t f = &ONE;
    FCLAW2D_FORT_INTEGRATE_EXACT(&mx,&my,&mbc,&dx,&dy,
                                 &xlower, &ylower, &blockno,
                                 area, &f, favg, &compute_avg,
                                 &ghost_only);

    FCLAW_FREE(favg);
}


void fclaw2d_metric_setup_mesh(fclaw2d_global_t *glob,
                               fclaw2d_patch_t *this_patch,
                               int blockno,
                               int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *area;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,
                                  &xp,&yp,&zp,&xd,&yd,&zd,&area);

    /* Compute centers and corners of mesh cell */
    fclaw2d_vt()->fort_setup_mesh(&mx,&my,&mbc,&xlower,&ylower,&dx,&dy,&blockno,
                       xp,yp,zp,xd,yd,zd);

}

void fclaw2d_metric_compute_normals(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    int blockno,
                                    int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *xnormals, *ynormals;
    double *xtangents, *ytangents;
    double *edgelengths;
    double *surfnormals, *curvature;
    double *area;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,
                                  &xp,&yp,&zp,&xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_data2(glob,this_patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,&edgelengths,
                                   &curvature);

    /* The user could set these to NULL to avoid doing these computations ... */

    if (fclaw2d_vt()->fort_compute_normals != NULL)
    {
        fclaw2d_vt()->fort_compute_normals(&mx,&my,&mbc,xp,yp,zp,xd,yd,zd,
                                xnormals,ynormals);
    }

    if (fclaw2d_vt()->fort_compute_tangents != NULL)
    {
        fclaw2d_vt()->fort_compute_tangents(&mx,&my,&mbc,xd,yd,zd,xtangents,ytangents,
                                 edgelengths);
    }

    if (fclaw2d_vt()->fort_compute_surf_normals != NULL)
    {
        fclaw2d_vt()->fort_compute_surf_normals(&mx,&my,&mbc,xnormals,ynormals,edgelengths,
                                     curvature, surfnormals, area);
    }

}
