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

#include "metric_user.h"
#include <fclaw2d_farraybox.hpp>

static fclaw2d_vtable_t fclaw2d_vt;

/* Want to write this so clawpack isn't needed */
static fc2d_clawpack46_vtable_t classic_claw46;

void metric_link_patch(fclaw_domain_t *domain)
{
    fclaw2d_init_vtable(&fclaw2d_vt);
    const user_options_t* user = metric_user_get_options(domain);

    fclaw2d_vt.problem_setup = &metric_problem_setup;

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt,&classic_claw46);

        fclaw2d_vt.patch_initialize = &metric_patch_initialize;
        fclaw2d_vt.metric_compute_area = &fclaw2d_metric_compute_area_exact;

        fclaw2d_vt.run_user_diagnostics = &metric_diagnostics;
        fclaw2d_vt.fort_compute_patch_error  = &compute_error;

        fc2d_clawpack46_set_vtable(classic_claw46);
    }

    fclaw2d_set_vtable(domain,&fclaw2d_vt);

}

void metric_problem_setup(fclaw_domain_t* domain)
{
    const user_options_t* user = metric_user_get_options(domain);

    /* Any general problem set up here */
    METRIC_SETPROB(&user->beta);  /* Set value of pi */
}

void metric_patch_initialize(fclaw_domain_t *domain,
                             fclaw_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *xnormals, *ynormals, *xtangents, *ytangents;
    double *surfnormals, *edgelengths;
    double *q, *area, *curvature;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    fclaw2d_clawpatch_metric_data2(domain, this_patch,
                                   &xnormals, &ynormals,
                                   &xtangents, &ytangents,
                                   &surfnormals, &edgelengths,
                                   &curvature);

    area =  fclaw2d_clawpatch_get_area(domain, this_patch);

#if 0
    error = fclaw2d_clawpatch_get_error(domain, this_patch);
    blockno = this_block_idx;
    compute_error(meqn,mbc,mx,my,blockno,xlower,ylower,dx,dy,
                  curvature,error);
#endif


    initialize(mx,my,meqn,mbc,xlower,ylower,dx,dy,q,
               curvature,area);
}
