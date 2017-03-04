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

#include "filament_user.h"

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>

#if 0
static fclaw2d_vtable_t fclaw2d_vt;

static fc2d_clawpack46_vtable_t classic_claw46;
static fc2d_clawpack5_vtable_t classic_claw5;
#endif 

void filament_link_solvers(fclaw2d_domain_t *domain)
{
    const user_options_t* user = filament_user_get_options(domain);
    const amr_options_t* gparms = fclaw2d_forestclaw_get_options(domain);

    // fclaw2d_init_vtable();

    if (user->claw_version == 4)
    {
        // fc2d_clawpack46_set_vtable_defaults(&fclaw2d_vt, &classic_claw46);
        // fc2d_clawpack46_set_vtable_defaults();
        
        fc2d_clawpack46_vt()->setprob   = &SETPROB;
        fc2d_clawpack46_vt()->qinit     = &CLAWPACK46_QINIT;
        if (gparms->manifold)
        {
            fclaw2d_patch_vt()->patch_setup   = &filament_patch_setup_manifold;

            fc2d_clawpack46_vt()->rpn2      = &CLAWPACK46_RPN2ADV_MANIFOLD;
            fc2d_clawpack46_vt()->rpt2      = &CLAWPACK46_RPT2ADV_MANIFOLD;

            if (user->example == 2)
            {
                /* Avoid tagging block corners in 5 patch example*/
                fclaw2d_clawpatch_vt()->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
                fclaw2d_clawpatch_vt()->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;
            }
        }
        else
        {
            fc2d_clawpack46_vt()->setaux    = &CLAWPACK46_SETAUX;  /* Used in non-manifold case */
            fc2d_clawpack46_vt()->rpn2      = &CLAWPACK46_RPN2ADV;
            fc2d_clawpack46_vt()->rpt2      = &CLAWPACK46_RPT2ADV;
        }

        // fc2d_clawpack46_set_vtable(classic_claw46);
    }
    else if (user->claw_version == 5)
    {
        // fc2d_clawpack5_set_vtable_defaults(&fclaw2d_vt, &classic_claw5);
        // fc2d_clawpack5_set_vtable_defaults();

        fc2d_clawpack5_vt()->setprob   = &SETPROB;
        fc2d_clawpack5_vt()->qinit     = &CLAWPACK5_QINIT;
        if (gparms->manifold)
        {
            fclaw2d_patch_vt()->patch_setup  = &filament_patch_setup_manifold;
            fc2d_clawpack5_vt()->rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
            fc2d_clawpack5_vt()->rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;

            if (user->example == 2)
            {
                /* Avoid tagging block corners in 5 patch example*/
                fclaw2d_clawpatch_vt()->fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;
                fclaw2d_clawpatch_vt()->fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;
            }
        }
        else
        {
            fc2d_clawpack5_vt()->setaux    = &CLAWPACK5_SETAUX;   /* Used in non-manifold case */
            fc2d_clawpack5_vt()->rpn2      = &CLAWPACK5_RPN2ADV;
            fc2d_clawpack5_vt()->rpt2      = &CLAWPACK5_RPT2ADV;
        }

        // fc2d_clawpack5_set_vtable(classic_claw5);
    }
    fclaw2d_set_vtable();
}

void filament_patch_setup_manifold(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *this_patch,
                                   int this_block_idx,
                                   int this_patch_idx)
{
    const user_options_t* user = filament_user_get_options(domain);

    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    if (user->claw_version == 4)
    {
        fc2d_clawpack46_define_auxarray(domain,this_patch);
        fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);
        USER46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                               &dx,&dy,&maux,aux,&this_block_idx,
                               xd,yd,zd,area);
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_define_auxarray(domain,this_patch);
        fc2d_clawpack5_aux_data(domain,this_patch,&aux,&maux);
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                              &dx,&dy,&maux,aux,&this_block_idx,
                              xd,yd,zd,area);
    }
}
