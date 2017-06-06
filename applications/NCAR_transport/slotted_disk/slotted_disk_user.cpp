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

#include "slotted_disk_user.h"

#include <fclaw2d_include_all.h>

#include "../all/transport_user.h"
#include "../all/transport_options.h"

void slotted_disk_link_solvers(fclaw2d_global_t *glob)
{
    transport_link_solvers(glob);

    /* Custom setprob */
    fclaw2d_vtable_t *vt = fclaw2d_vt();
    vt->problem_setup    = &slotted_disk_problem_setup;  /* Version-independent */
}


void slotted_disk_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = transport_get_options(glob);
    const fclaw_options_t* gparms = fclaw2d_get_options(glob);
    SLOTTED_DISK_SETPROB(&user->kappa, &gparms->tfinal);
}


#if 0
void slotted_disk_patch_setup(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    const user_options_t* user = slotted_disk_user_get_options(domain);

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        return;
    }

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_define_auxarray(domain,this_patch);
        fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

        USER46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                               &maux,aux,&this_block_idx,xd,yd,zd,area);
    }
    else if(user->claw_version == 5)
    {
        fc2d_clawpack5_define_auxarray(domain,this_patch);
        fc2d_clawpack5_aux_data(domain,this_patch,&aux,&maux);

        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                              &maux,aux,&this_block_idx,xd,yd,zd,area);
    }
}

void slotted_disk_b4step2(fclaw2d_domain_t *domain,
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
    const user_options_t* user = slotted_disk_user_get_options(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,&xd,&yd,&zd,&area);

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);
        USER46_B4STEP2_MANIFOLD(&mx,&my,&mbc,&dx,&dy,&t,&maux,aux,&this_block_idx,xd,yd,zd);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_aux_data(domain,this_patch,&aux,&maux);
        USER5_B4STEP2_MANIFOLD(&mx,&my,&mbc,&dx,&dy,&t,&maux,aux,&this_block_idx,xd,yd,zd);
    }
}

double slotted_disk_update(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt)
{
    const user_options_t* user = slotted_disk_user_get_options(domain);

    slotted_disk_b4step2(domain,this_patch,this_block_idx,
                         this_patch_idx,t,dt);

    double maxcfl;
    if (user->claw_version == 4)
    {
        maxcfl = fc2d_clawpack46_step2(domain,this_patch,this_block_idx,
                                       this_patch_idx,t,dt);
    }
    else
    {
        maxcfl = fc2d_clawpack5_step2(domain,this_patch,this_block_idx,
                                      this_patch_idx,t,dt);
    }

    return maxcfl;
}
#endif
