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

#include "amr_forestclaw.h"
#include "forestclaw2d.h"
#include "fc2d_clawpack46.H"
#include "swirl_user.H"
#include "fclaw2d_vtable.h"
#include "fclaw2d_output_ascii.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* Most of these do not use the "classic" signature. */

static fclaw2d_vtable_t vt;
static fc2d_clawpack46_vtable_t classic_user;

void swirl_link_solvers(fclaw2d_domain_t *domain)
{
    vt.problem_setup            = &fc2d_clawpack46_setprob;

    vt.patch_setup              = &swirl_patch_setup;
    vt.patch_initialize         = &fc2d_clawpack46_qinit;
    vt.patch_physical_bc        = &fc2d_clawpack46_bc2;
    vt.patch_single_step_update = &swirl_patch_single_step_update;

    vt.patch_tag4refinement     = &swirl_patch_tag4refinement;
    vt.patch_tag4coarsening     = &swirl_patch_tag4coarsening;

    vt.write_header             = &fclaw2d_output_header_ascii;
    vt.write_tfile              = &FCLAW2D_OUTPUT_WRITE_TFILE;

    vt.patch_write_file         = &fclaw2d_output_patch_ascii;
    vt.patch_write_qfile        = &FCLAW2D_OUTPUT_WRITE_QFILE;

    fclaw2d_set_vtable(domain,&vt);

    /* Needed for the clawpack46 package */
    classic_user.setprob = &SETPROB;
    classic_user.qinit = &QINIT;
    classic_user.rpn2 = &RPN2;
    classic_user.rpt2 = &RPT2;

    fc2d_clawpack46_set_vtable(&classic_user);
}

void swirl_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fc2d_clawpack46_define_auxarray2(domain,this_patch);

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    SETAUX_MANIFOLD(mbc,mx,my,xlower,ylower,dx,dy,maux,
                    aux,this_block_idx,xd,yd,zd,area);
}

void swirl_patch_b4step2(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    /* Update the velocity field */
    B4STEP2_MANIFOLD(mbc,mx,my,dx,dy,this_block_idx,
                     xd,yd,zd,t, dt,maux,aux);
}


double swirl_patch_single_step_update(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double dt)
{
    swirl_patch_b4step2(domain,this_patch,this_block_idx,
                        this_patch_idx,t,dt);

    double maxcfl = fc2d_clawpack46_step2(domain,this_patch,this_block_idx,
                                          this_patch_idx,t,dt);
    return maxcfl;
}


fclaw_bool swirl_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno, int this_patch_idx,
                                      int initflag)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    tag_patch = 0;
    swirl_tag4refinement_(mx,my,mbc,meqn,blockno,xlower,ylower,
                          dx,dy,q,initflag,tag_patch);
    return tag_patch == 1;
}

fclaw_bool swirl_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qcoarse;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&qcoarse,&meqn);

    tag_patch = 1;
    swirl_tag4coarsening_(mx,my,mbc,meqn,xlower,ylower,dx,dy,qcoarse,tag_patch);
    return tag_patch == 0;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
