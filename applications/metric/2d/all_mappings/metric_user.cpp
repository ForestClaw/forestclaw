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

#include "amr_forestclaw.H"
#include "ClawPatch.H"
#include "fclaw2d_map_query.h"
#include "fclaw2d_clawpatch.h"
#include "fclaw2d_physical_bc.h"
#include "fclaw2d_vtable.h"

#include "metric_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fclaw2d_vtable_t vt;

void metric_link_patch(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);

    vt.problem_setup = &metric_problem_setup;

    vt.patch_initialize = &metric_patch_initialize;
    vt.patch_physical_bc = &fclaw2d_physical_bc_default;  /* Doesn't do anything */

    vt.run_diagnostics = &metric_diagnostics;

    fclaw2d_set_vtable(domain,&vt);

}

void metric_problem_setup(fclaw2d_domain_t* domain)
{
    /* Any general problem set up here */
    SETPROB();  /* Set value of pi */
}

void metric_patch_initialize(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q, *area, *curvature, *error_ptr;
    ClawPatch *cp;
    int blockno;

    cp = get_clawpatch(this_patch);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    area = cp->area();
    curvature = cp->curvature();

    /* Create an array with same dimensions as q, and one field */
    FArrayBox error;
    error.define(cp->dataBox(),1);
    error_ptr = error.dataPtr();

    fclaw2d_map_context_t* cont = get_map_context(domain);
    blockno = this_block_idx;
    compute_error(meqn,mbc,mx,my,&cont,blockno,xlower,ylower,dx,dy,
                  curvature,error_ptr);
    initialize(mx,my,meqn,mbc,xlower,ylower,dx,dy,q,
               error_ptr,curvature,area);
}



#ifdef __cplusplus
#if 0
{
#endif
}
#endif
