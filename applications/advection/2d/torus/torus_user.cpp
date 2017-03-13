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

#include "torus_user.h"

#include <fclaw2d_forestclaw.h>
#include "fclaw2d_clawpatch.h"

void torus_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t *user =  torus_user_get_options(domain);

    fclaw2d_vt()->problem_setup = &torus_problem_setup;
    fclaw2d_patch_vt()->setup = &torus_patch_setup;

    if (user->claw_version == 4)
    {
        fclaw2d_clawpatch_vt()->fort_compute_patch_error = &TORUS46_COMPUTE_ERROR;

        fc2d_clawpack46_vt()->qinit     = &CLAWPACK46_QINIT;
        fc2d_clawpack46_vt()->setaux    = &TORUS46_SETAUX;  /* Not really a mapped setaux */
        fc2d_clawpack46_vt()->rpn2      = &CLAWPACK46_RPN2ADV_MANIFOLD;
        fc2d_clawpack46_vt()->rpt2      = &CLAWPACK46_RPT2ADV_MANIFOLD;

        if (user->example == 1)
        {
            /* Accuracy problem : Used divided differences for tagging */
            fclaw2d_clawpatch_vt()->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
            fclaw2d_clawpatch_vt()->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

            /* Include error in output files */
            fclaw2d_clawpatch_vt()->fort_write_header   = &TORUS_FORT_WRITE_HEADER;
            fclaw2d_patch_vt()->write_file        = &torus_output_write_file;
        }

    }
    else if (user->claw_version == 5)
    {
        fclaw2d_clawpatch_vt()->fort_compute_patch_error = &TORUS5_COMPUTE_ERROR;

        fclaw2d_patch_vt()->setup      = &torus_patch_setup;

        fc2d_clawpack5_vt()->qinit     = &CLAWPACK5_QINIT;
        fc2d_clawpack5_vt()->setaux    = &TORUS5_SETAUX;
        fc2d_clawpack5_vt()->rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        fc2d_clawpack5_vt()->rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;

        if (user->example == 1)
        {
            /* Accuracy problem : Used divided differences for tagging */
            fclaw2d_clawpatch_vt()->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
            fclaw2d_clawpatch_vt()->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

            /* Write out error */
            fclaw2d_clawpatch_vt()->fort_write_header      = &TORUS_FORT_WRITE_HEADER;
            fclaw2d_patch_vt()->write_file    = &torus_output_write_file;
        }
    }
}

void torus_problem_setup(fclaw2d_domain_t *domain)
{
    const user_options_t* user = torus_user_get_options(domain);
    TORUS_SETPROB(&user->example,&user->alpha);
}

void torus_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    const user_options_t* user = torus_user_get_options(domain);

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_setaux(domain,this_patch,this_block_idx,this_patch_idx);
        fc2d_clawpack46_set_capacity(domain,this_patch,this_block_idx,this_patch_idx);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_setaux(domain,this_patch,this_block_idx,this_patch_idx);
        fc2d_clawpack5_set_capacity(domain,this_patch,this_block_idx,this_patch_idx);
    }
}

void torus_output_write_file(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx, int this_patch_idx,
                             int iframe, int patch_num,int level)
{
    /* This new wrapper is needed because we are passing both q
       and the error into the FORT file.  */
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy,t;
    double *q, *error;
    char matname1[11];

    const user_options_t *user = torus_user_get_options(domain);

    t = fclaw2d_domain_get_time(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    error = fclaw2d_clawpatch_get_error(domain,this_patch);

    sprintf(matname1,"fort.q%04d",iframe);

    /* Here, we pass in q and the error, so need special headers and files */
    if (user->claw_version == 4)
    {
        TORUS46_FORT_WRITE_FILE(matname1, &mx,&my,&meqn,&mbc,&xlower,&ylower,
                                &dx,&dy,q,error,&t,
                                &patch_num,&level,&this_block_idx,
                                &domain->mpirank);
    }
    else if (user->claw_version == 5)
    {
        TORUS5_FORT_WRITE_FILE(matname1, &mx,&my,&meqn,&mbc,&xlower,&ylower,
                               &dx,&dy,q,error,&t,
                               &patch_num,&level,&this_block_idx,
                               &domain->mpirank);
    }
}
