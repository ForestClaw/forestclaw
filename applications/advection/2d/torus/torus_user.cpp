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

#include "amr_includes.H"
#include "fc2d_clawpack46.H"
#include "torus_user.H"

#include <fclaw2d_vtable.h>
#include <fclaw_register.h>
#include <fclaw2d_output.h>
#include <fclaw2d_output_ascii.h>
#include <fclaw2d_regrid_default.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fc2d_clawpack46_vtable_t classic_claw;

static fclaw2d_vtable_t vt;

void torus_link_solvers(fclaw2d_domain_t *domain)
{
    const amr_options_t *gparms;
    int m;

    fclaw2d_init_vtable(&vt);

    gparms = fclaw2d_forestclaw_get_options(domain);
    m = gparms->manifold;

    vt.problem_setup            = &fc2d_clawpack46_setprob;

    vt.patch_setup = m ? &torus_patch_manifold_setup : &fc2d_clawpack46_setaux;

    vt.patch_initialize         = &fc2d_clawpack46_qinit;
    vt.patch_physical_bc        = &fc2d_clawpack46_bc2;     /* Needed for lat-long grid */
    vt.patch_single_step_update = &fc2d_clawpack46_update;  /* Includes b4step2 and src2 */

    vt.patch_tag4refinement     = &fclaw2d_patch_tag4refinement;
    vt.fort_tag4refinement      = &FCLAW2D_FORT_TAG4REFINEMENT;

    vt.patch_tag4coarsening     = &fclaw2d_patch_tag4coarsening;
    vt.fort_tag4coarsening      = &FCLAW2D_FORT_TAG4COARSENING;

    vt.write_header             = &fclaw2d_output_header_ascii;
    vt.fort_write_header        = &FCLAW2D_FORT_WRITE_HEADER;

    vt.patch_write_file         = &fclaw2d_output_patch_ascii;
    vt.fort_write_file          = &FCLAW2D_FORT_WRITE_FILE;

    fclaw2d_set_vtable(domain,&vt);

    /* Needed for the clawpack46 package */
    classic_claw.qinit = &QINIT;
    classic_claw.setaux = &SETAUX;    /* Called by fc2d_clawpack46_setaux, above */
    classic_claw.setprob = &SETPROB;
    classic_claw.rpn2 = &RPN2;
    classic_claw.rpt2 = &RPT2;

    fc2d_clawpack46_set_vtable(&classic_claw);

}

void torus_patch_manifold_setup(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *aux;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fc2d_clawpack46_define_auxarray2(domain,this_patch);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    SETAUX_MANIFOLD(&mx,&my,&mbc,&xlower,&ylower,&dx,&dy,
                    &this_block_idx, &maux, aux, area);
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
