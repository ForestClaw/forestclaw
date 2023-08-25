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
#include <fclaw_forestclaw.h>

#include <fc2d_clawpack46.h>


static fc2d_clawpack46_vtable_t classic_claw;

static fclaw2d_vtable_t vt;

void sphere_link_solvers(fclaw_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);
    fc2d_clawpack46_init_vtable(&classic_claw);

    vt.problem_setup = &fc2d_clawpack46_setprob;
    classic_claw.setprob = &SETPROB;

    vt.patch_setup = &sphere_patch_manifold_setup;
    /* classic_claw.setaux = &SETAUX_SPHERE; */

    vt.patch_initialize = &fc2d_clawpack46_qinit;
    classic_claw.qinit = &QINIT;

    vt.patch_physical_bc = &fc2d_clawpack46_bc2;     /* Needed for lat-long grid */

    vt.metric_compute_area = &fclaw2d_metric_compute_area;

    vt.patch_single_step_update = &fc2d_clawpack46_update;  /* Includes b4step2 and src2 */
    classic_claw.b4step2 = &B4STEP2;
    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;


    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);

}

void sphere_patch_manifold_setup(fclaw_domain_t *domain,
                                fclaw_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux;
    double *xnormals, *ynormals,*xtangents;
    double *ytangents,*surfnormals;
    double *edgelengths, *curvature;

    fc2d_clawpack46_define_auxarray(domain,this_patch);

    fclaw_clawpatch_grid_data_2d(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);
    fclaw2d_clawpatch_metric_data2(domain,this_patch,&xnormals,&ynormals,
                                   &xtangents,&ytangents,&surfnormals,
                                   &edgelengths,&curvature);

    SETAUX_SPHERE(&mx,&my,&mbc,&xlower,&ylower,
                  &dx,&dy,&maux,aux,xnormals,ynormals,
                  xtangents,ytangents,surfnormals);

    fc2d_clawpack46_set_capacity(domain,this_patch,this_block_idx,this_patch_idx);
}
