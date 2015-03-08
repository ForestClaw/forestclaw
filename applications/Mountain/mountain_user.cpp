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

#include <amr_includes.H>
#include <fc2d_clawpack46.H>
#include <fclaw2d_vtable.h>
#include <fclaw2d_regrid.h>
#include "mountain_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fc2d_clawpack46_vtable_t classic_claw;
static fclaw2d_vtable_t vt;

void mountain_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);
    fc2d_clawpack46_init_vtable(&classic_claw);

    vt.problem_setup = &fc2d_clawpack46_setprob;
    classic_claw.setprob = &SETPROB;

    vt.patch_setup = &mountain_patch_setup;
    classic_claw.setaux = &SETAUX;  /* Called by fc2d_clawpack46_setaux */

    vt.patch_initialize = &fc2d_clawpack46_qinit;
    classic_claw.qinit = &QINIT;

    vt.patch_physical_bc = &fc2d_clawpack46_bc2;   /* This doesn't do anything */

    /* Use this to refine near the terrain.  The default is to refine near
       advection front */
    vt.fort_tag4refinement = &TAG4REFINEMENT; /* Customized refinement criteria */
    vt.fort_tag4coarsening = &TAG4COARSENING; /* Customized refinement criteria */

    vt.patch_single_step_update = &fc2d_clawpack46_update;
    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;

    fclaw2d_set_vtable(domain,&vt);
    fc2d_clawpack46_set_vtable(&classic_claw);
}

void mountain_patch_setup(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_blockno,
                          int this_patchno)
{
    fc2d_clawpack46_setaux(domain,this_patch,this_blockno,
                           this_patchno);
    fc2d_clawpack46_set_capacity(domain,this_patch,this_blockno,
                                 this_patchno);
}

double mountain_patch_update(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    /* Don't do anything */
    return gparms->desired_cfl;
}




#ifdef __cplusplus
#if 0
{
#endif
}
#endif
