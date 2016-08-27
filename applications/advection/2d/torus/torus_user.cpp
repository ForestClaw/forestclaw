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

#include <fclaw2d_forestclaw.h>
#include "torus_common.h"
#include "torus_user.h"
#include "fclaw2d_clawpatch.h"

#include <fc2d_clawpack46.h>

static fc2d_clawpack46_vtable_t classic_claw;

static fclaw2d_vtable_t vt;

void torus_link_solvers(fclaw2d_domain_t *domain)
{

    fclaw_app_t* app = fclaw2d_domain_get_app(domain);
    user_options_t* user = (user_options_t*) fclaw_app_get_user(app);
    int example = user->example;
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_init_vtable(&vt);
    fc2d_clawpack46_init_vtable(&vt, &classic_claw);

    vt.problem_setup            = &torus_patch_setup;

    if (gparms->manifold)
    {
        vt.patch_setup = &torus_patch_manifold_setup;
        classic_claw.setaux = &SETAUX_MANIFOLD;
    }
    else
    {
        vt.patch_setup = &fc2d_clawpack46_setaux;
        classic_claw.setaux = &SETAUX;
    }


    vt.patch_initialize         = &fc2d_clawpack46_qinit;
    classic_claw.qinit = &QINIT;

    vt.patch_physical_bc        = &fc2d_clawpack46_bc2;     /* Needed for lat-long grid */

    vt.patch_single_step_update = &fc2d_clawpack46_update;  /* Includes b4step2 and src2 */
    if (gparms->manifold)
    {
        classic_claw.b4step2 = &B4STEP2_MANIFOLD;
    }
    else
    {
        /* classic_claw.b4step2 = &B4STEP2; */
    }

    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;

    vt.fort_compute_patch_error = &TORUS_COMPUTE_ERROR;

    if (example == 6)
    {
        vt.fort_tag4refinement = &TAG4REFINEMENT;
        vt.fort_tag4coarsening = &TAG4COARSENING;

        vt.write_header      = &torus_output_write_header;
        vt.patch_write_file  = &torus_output_write_file;
    }

    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);

}

void torus_patch_setup(fclaw2d_domain_t *domain)
{
    fclaw_app_t* app = fclaw2d_domain_get_app(domain);
    user_options_t* user = (user_options_t*) fclaw_app_get_user(app);
    int example = user->example;

    SETPROB_TORUS(&example);
}

void torus_patch_manifold_setup(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{

    fc2d_clawpack46_setaux(domain,this_patch,this_block_idx,this_patch_idx);
    fc2d_clawpack46_set_capacity(domain,this_patch,this_block_idx,this_patch_idx);
}


void torus_output_write_header(fclaw2d_domain_t* domain,
                               int iframe)
{
    const amr_options_t *amropt;
    fclaw2d_vtable_t vt;
    int meqn,ngrids;
    double time;
    char matname1[10];
    char matname2[10];

    amropt = fclaw2d_forestclaw_get_options(domain);

    time = fclaw2d_domain_get_time(domain);
    ngrids = fclaw2d_domain_get_num_patches(domain);

    meqn = amropt->meqn;

    sprintf(matname1,"fort.q%04d",iframe);
    sprintf(matname2,"fort.t%04d",iframe);

    vt = fclaw2d_get_vtable(domain);
    TORUS_FORT_WRITE_HEADER(matname1,matname2,&time,&meqn,&ngrids);
}


void torus_output_write_file(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx, int this_patch_idx,
                             int iframe, int patch_num,int level)
{
    fclaw2d_vtable_t vt;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy,t;
    double *q, *error;
    char matname1[10];
    vt = fclaw2d_get_vtable(domain);

    t = fclaw2d_domain_get_time(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    error = fclaw2d_clawpatch_get_error(domain,this_patch);

    sprintf(matname1,"fort.q%04d",iframe);
    TORUS_FORT_WRITE_FILE(matname1, &mx,&my,&meqn,&mbc,&xlower,&ylower,
                          &dx,&dy,q,error,&t,
                          &patch_num,&level,&this_block_idx,
                          &domain->mpirank);
}
