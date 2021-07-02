/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#include "hemisphere_user.h"

#include "../all/advection_user.h"

void hemisphere_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    patch_vt->setup      = &hemisphere_patch_setup;    

    const user_options_t* user = hemisphere_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();

        claw46_vt->fort_setprob   = &SETPROB;
        claw46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV_MANIFOLD;
        claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV_MANIFOLD;

        if (user->example == 0)
        {
            /* Avoid tagging block corners in 5 patch example */
            fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
            clawpatch_vt->fort_tag4refinement = &CLAWPATCH46_TAG4REFINEMENT;
            clawpatch_vt->fort_tag4coarsening = &CLAWPATCH46_TAG4COARSENING;
        }

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt();

        claw5_vt->fort_setprob   = &SETPROB;
        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        claw5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV_MANIFOLD;
        claw5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV_MANIFOLD;
        
        if (user->example == 0)
        {
            /* Avoid tagging block corners in 5 patch example */
            fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
            clawpatch_vt->fort_tag4refinement = &CLAWPATCH5_TAG4REFINEMENT;
            clawpatch_vt->fort_tag4coarsening = &CLAWPATCH5_TAG4COARSENING;
        }
    }
}


void hemisphere_patch_setup(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *patch,
                            int blockno,
                            int patchno)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    int maux;
    double *aux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);

    const user_options_t* user = hemisphere_get_options(glob);
    if (user->claw_version == 4)
    {
        USER46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                               &maux,aux,&blockno,xd,yd,zd,area);
    }
    else if (user->claw_version == 5)
    {
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                              &maux,aux,&blockno,xd,yd,zd,area);
    }

}
