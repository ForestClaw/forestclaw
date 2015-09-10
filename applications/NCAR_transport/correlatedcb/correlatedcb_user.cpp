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

#include "correlatedcb_user.h"

#include "fclaw2d_forestclaw.h"
#include "fclaw2d_clawpatch.h"
#include "fc2d_clawpack46.h"

static fclaw2d_vtable_t vt;
static fc2d_clawpack46_vtable_t classic_claw;

void correlatedcb_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);
    fc2d_clawpack46_init_vtable(&classic_claw);

    vt.problem_setup             = &correlatedcb_problem_setup;

    vt.patch_setup               = &correlatedcb_patch_setup;

    vt.patch_initialize          = &correlatedcb_qinit;
    vt.patch_physical_bc         = &fc2d_clawpack46_bc2;

    vt.patch_single_step_update  = &correlatedcb_update;
    classic_claw.rpn2            = &RPN2;
    classic_claw.rpt2            = &RPT2;

    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);
}


void correlatedcb_problem_setup(fclaw2d_domain_t* domain)
{
    const user_options_t *user;
    user = (user_options_t*) fclaw2d_domain_get_user_options(domain);
    SETPROB_TRANSPORT(&user->vflag,&user->init_choice);
}

void correlatedcb_patch_setup(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx)
{
    int mx, my, mbc, maux;
    double xlower,ylower, dx,dy;
    double *xp,*yp,*zp,*xd,*yd,*zd;
    double *aux, *area;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_clawpack46_define_auxarray(domain,this_patch);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    SETAUX_TRANSPORT(&mx,&my,&mbc,&xlower,&ylower,&dx,&dy,
                      &maux,aux,&this_block_idx,xd,yd,zd,area);
}


void correlatedcb_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    int mx, my, mbc, meqn, maux;
    double xlower,ylower, dx,dy;
    double *xp,*yp,*zp,*xd,*yd,*zd;
    double *aux, *area;
    double *q;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,&xd,&yd,&zd,&area);
    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    QINIT_TRANSPORT(&mx,&my,&meqn,&mbc,&xlower,&ylower,&dx,&dy,q,&maux,aux,
                     &this_block_idx,xp,yp,zp);
}


void correlatedcb_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    int mx, my, mbc, maux;
    double xlower,ylower, dx,dy;
    double *xp,*yp,*zp,*xd,*yd,*zd;
    double *aux, *area;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,&xd,&yd,&zd,&area);

    B4STEP2_TRANSPORT(&mx,&my,&mbc,&dx,&dy,&t,&maux,aux,&this_block_idx,xd,yd,zd);
}

double correlatedcb_update(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt)
{
    correlatedcb_b4step2(domain,this_patch,this_block_idx,
                         this_patch_idx,t,dt);

    double maxcfl = fc2d_clawpack46_step2(domain,this_patch,this_block_idx,
                                           this_patch_idx,t,dt);

    return maxcfl;
}
