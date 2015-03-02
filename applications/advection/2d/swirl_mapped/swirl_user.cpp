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

#include "forestclaw2d.h"
#include "fc2d_clawpack46.H"
#include "fclaw2d_vtable.h"
#include "fclaw2d_output_ascii.h"
#include "fclaw2d_regrid_default.h"

#include "swirl_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* Most of these do not use the "classic_claw" signature. */

static fclaw2d_vtable_t vt;
static fc2d_clawpack46_vtable_t classic_claw;

void swirl_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);

    vt.problem_setup            = &fc2d_clawpack46_setprob;

    vt.patch_setup              = &swirl_patch_setup;
    vt.patch_initialize         = &fc2d_clawpack46_qinit;
    vt.patch_physical_bc        = &fc2d_clawpack46_bc2;
    vt.patch_single_step_update = &swirl_patch_single_step_update;

    vt.patch_tag4refinement     = &fclaw2d_patch_tag4refinement;
    vt.fort_tag4refinement      = &FCLAW2D_FORT_TAG4REFINEMENT;

    vt.patch_tag4coarsening     = &fclaw2d_patch_tag4coarsening;
    vt.fort_tag4coarsening      = &FCLAW2D_FORT_TAG4COARSENING;

    vt.write_header             = &fclaw2d_output_header_ascii;
    vt.fort_write_header        = &FCLAW2D_FORT_WRITE_HEADER;

    vt.patch_write_file         = &fclaw2d_output_patch_ascii;
    vt.fort_write_file          = &FCLAW2D_FORT_WRITE_FILE;

    fclaw2d_set_vtable(domain,&vt);

    /* Needed for the clawpack46 package */
    classic_claw.setprob = &SETPROB;
    classic_claw.qinit = &QINIT;
    classic_claw.bc2 = NULL;      /* Use clawpack46 library routine */
    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;

    fc2d_clawpack46_set_vtable(&classic_claw);
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



#ifdef __cplusplus
#if 0
{
#endif
}
#endif
