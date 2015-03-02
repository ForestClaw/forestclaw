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

#include <fclaw2d_vtable.h>
#include <fclaw_register.h>
#include <fclaw2d_output_ascii.h>
#include <fclaw2d_regrid_default.h>

#include "no_solver_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fclaw2d_vtable_t vt;


void no_solver_linker(fclaw2d_domain_t* domain)
{
    fclaw2d_init_vtable(&vt);

    vt.patch_initialize         = &no_solver_patch_initialize;
    vt.patch_physical_bc        = &fclaw2d_physbc_default;  /* No BCs are imposed */
    vt.patch_single_step_update = &no_solver_update;

    fclaw2d_set_vtable(domain,&vt);
}

void no_solver_patch_initialize(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    int blockno = this_block_idx;
    INITIALIZE(&mx,&my,&meqn,&mbc,&blockno,&xlower,&ylower,&dx,&dy,q);
}

double no_solver_update(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt)
{
    const amr_options_t *gparms;
    gparms = fclaw2d_forestclaw_get_options(domain);

    /* This is needed to avoid a floating point errror */
    ClawPatch *cp = get_clawpatch(this_patch);
    cp->save_current_step();  // Save for time interpolation

    return gparms->desired_cfl;
}



#ifdef __cplusplus
#if 0
{
#endif
}
#endif
