/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "no_solver_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_output_ascii.h>
#include <fclaw2d_physical_bc.h>

#include "gem3d_output_mesh.h"

static
void no_solver_output(fclaw2d_global_t* glob,int iframe);

void no_solver_link_solvers(fclaw2d_global_t* global)
{
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    fclaw2d_vtable_t* vt = fclaw2d_vt(glob);

    patch_vt->single_step_update   = &no_solver_update;
    patch_vt->initialize           = &no_solver_patch_initialize;
    patch_vt->physical_bc          = &fclaw2d_physical_bc_default;  /* do nothing */
    vt->output_frame               = &no_solver_output;

}

void no_solver_patch_initialize(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    int blockno = this_block_idx;
    INITIALIZE(&mx,&my,&meqn,&mbc,&blockno,&xlower,&ylower,&dx,&dy,q);
}

double no_solver_update(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt)
{
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    return fclaw_opt->desired_cfl;
}

static
void no_solver_output(fclaw2d_global_t* glob,int iframe)
{
    /* Create usual output */
    fclaw2d_clawpatch_output_ascii(glob,iframe);

    /* Write out mesh (for GEM3d) */
    gem3d_output_mesh(glob,iframe);
}

