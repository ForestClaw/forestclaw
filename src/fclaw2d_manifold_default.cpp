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
#include <fclaw2d_manifold_default_fort.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void fclaw2d_manifold_compute_area(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t* this_patch,
                                   int blockno,
                                   int patchno)
{
    fclaw2d_vtable_t vt;
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    int level, maxlevel,refratio;

    const amr_options_t* gparms = get_domain_parms(domain);
    level = this_patch->level;
    maxlevel = gparms->maxlevel;
    refratio = gparms->refratio;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *area = fclaw2d_clawpatch_get_area(domain,this_patch);

    vt = fclaw2d_get_vtable(domain);

    /* vt.fort_compute_area(...) */
    compute_area_(&mx, &my, &mbc, &dx, &dy, &xlower, &ylower,
                  &blockno, area, &level, &maxlevel, &refratio);
}


void fclaw2d_manifold_setup_mesh(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *this_patch,
                                 int blockno,
                                 int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *area;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,
                                  &xp,&yp,&zp,&xd,&yd,&zd,&area);

    /* Compute centers and corners of mesh cell */
    setup_mesh_(&mx,&my,&mbc,&xlower,&ylower,&dx,&dy,&blockno,
                xp,yp,zp,xd,yd,zd);

}

void fclaw2d_manifold_compute_normals(fclaw2d_domain_t *domain,
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

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,
                                  &xp,&yp,&zp,&xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_data2(domain,this_patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &edgelengths,
                                   &surfnormals,&curvature);


    /* vt.fort_compute_normals(...) */
    compute_normals_(&mx,&my,&mbc,xp,yp,zp,xd,yd,zd,
                     xnormals,ynormals);

    compute_tangents_(&mx,&my,&mbc,xd,yd,zd,xtangents,ytangents,edgelengths);

}

void fclaw2d_manifold_compute_curvature(fclaw2d_domain_t *domain,
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

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,
                                  &xp,&yp,&zp,&xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_data2(domain,this_patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &edgelengths,
                                   &surfnormals,&curvature);

    /* vt.fort_compute_curvature(...) */
    compute_surf_normals_(&mx,&my,&mbc,xnormals,ynormals,edgelengths,
                          curvature, surfnormals, area);


}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
