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

#include "radial_user.h"

static
void radial_problem_setup(fclaw_global_t* glob)
{
    user_options_t* user = radial_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24.16f   %s",user->rho,"\% rho\n");
        fprintf(f,  "%-24.16f   %s",user->bulk,"\% bulk\n");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
    fclaw2d_domain_barrier (glob->domain);  /* redundant?  */
 
    /* rho, bulk are inputs; cc and zz are outputs.  Results are
       stored in a common block */
    SETPROB();
}

static
void radial_patch_setup(fclaw_global_t *glob,
                        fclaw_patch_t *patch,
                        int blockno,
                        int patchno)
{
    if (fclaw2d_patch_is_ghost(patch))
    {
        /* Mapped info is needed only for an update */
        return;
    }

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xp,*yp,*zp;
    double *xd,*yd,*zd,*area;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *xnormals,*ynormals,*xtangents,*ytangents;
    double *surfnormals,*edgelengths,*curvature;
    fclaw2d_clawpatch_metric_data2(glob,patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,&edgelengths,
                                   &curvature);

    int maux;
    double* aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    const user_options_t* user = radial_get_options(glob);
    
    if (user->claw_version == 4)
        CLAWPACK46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                                   &dx,&dy,&maux,aux,
                                   xnormals,ynormals,edgelengths,area);   

    else if (user->claw_version == 5)    
        CLAWPACK5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                                  &dx,&dy,&maux,aux,
                                  xnormals,ynormals,edgelengths,area);

}

void radial_link_solvers(fclaw_global_t *glob)
{
    fclaw_vtable_t *vt = fclaw_vt(glob);

    vt->problem_setup = &radial_problem_setup;  /* Version-independent */

    const user_options_t* user = radial_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt(glob);

        claw46_vt->fort_qinit     = &CLAWPACK46_QINIT;

        if (user->example == 0)
        {
            claw46_vt->fort_rpn2      = &CLAWPACK46_RPN2;
            claw46_vt->fort_rpt2      = &CLAWPACK46_RPT2;
        }
        else if (user->example == 1 || user->example == 2)
        {
            fclaw2d_patch_vtable_t  *patch_vt = fclaw2d_patch_vt(glob);
            
            patch_vt->setup = &radial_patch_setup;

            claw46_vt->fort_rpn2  = &CLAWPACK46_RPN2_MANIFOLD;
            claw46_vt->fort_rpt2  = &CLAWPACK46_RPT2_MANIFOLD;
        }
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t    *claw5_vt = fc2d_clawpack5_vt(glob);

        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;

        if (user->example == 0)
        {
            claw5_vt->fort_rpn2 = &CLAWPACK5_RPN2;
            claw5_vt->fort_rpt2 = &CLAWPACK5_RPT2;
        }
        else if (user->example == 1 || user->example == 2)
        {
            fclaw2d_patch_vtable_t  *patch_vt = fclaw2d_patch_vt(glob);
            
            patch_vt->setup = &radial_patch_setup;

            claw5_vt->fort_rpn2  = &CLAWPACK5_RPN2_MANIFOLD;
            claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2_MANIFOLD;
        }
    }
}


