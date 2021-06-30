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

#include "radialdam_user.h"

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

#include "../rp/shallow_user_fort.h"

void radialdam_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();

    vt->problem_setup = &radialdam_problem_setup;  /* Version-independent */

    const user_options_t* user = radialdam_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();
        claw46_vt->fort_qinit     = &CLAWPACK46_QINIT;

        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        clawpatch_vt->fort_tag4refinement = &CLAWPATCH46_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &CLAWPATCH46_TAG4COARSENING;

        if (user->example == 0) 
        {
            claw46_vt->fort_rpn2 = &CLAWPACK46_RPN2;
            claw46_vt->fort_rpt2 = &CLAWPACK46_RPT2;
            claw46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE;
        }
        else if (user->example >= 1 && user->example <= 3)
        {
            fclaw2d_patch_vtable_t  *patch_vt = fclaw2d_patch_vt();
            patch_vt->setup = &radialdam_patch_setup;

            claw46_vt->fort_rpn2  = &CLAWPACK46_RPN2_MANIFOLD;
            claw46_vt->fort_rpt2  = &CLAWPACK46_RPT2_MANIFOLD;
            claw46_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;
        }

    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt();
        claw5_vt->fort_qinit = &CLAWPACK5_QINIT;

        /* Avoid tagging block corners in 5 patch example*/
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        clawpatch_vt->fort_tag4refinement = &CLAWPATCH5_TAG4REFINEMENT;
        clawpatch_vt->fort_tag4coarsening = &CLAWPATCH5_TAG4COARSENING;
            
        if (user->example == 0)
        {
            claw5_vt->fort_rpn2 = &CLAWPACK5_RPN2;
            claw5_vt->fort_rpt2 = &CLAWPACK5_RPT2;
            claw5_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE;
        }
        else if (user->example >= 1 && user->example <= 3)
        {
            fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
            patch_vt->setup = &radialdam_patch_setup;

            claw5_vt->fort_rpn2  = &CLAWPACK5_RPN2_MANIFOLD;
            claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2_MANIFOLD;
            claw5_vt->fort_rpn2_cons = &RPN2_CONS_UPDATE_MANIFOLD;
        }
    }
}


void radialdam_problem_setup(fclaw2d_global_t* glob)
{
    const user_options_t* user = radialdam_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",   user->example,"\% example\n");
        fprintf(f,  "%-24.16f   %s",user->g,"\% gravity\n");
        fprintf(f,  "%-24.16f   %s",user->x0,"\% x0\n");
        fprintf(f,  "%-24.16f   %s",user->y0,"\% y0\n");
        fprintf(f,  "%-24.16f   %s",user->r0,"\% r0\n");
        fprintf(f,  "%-24.16f   %s",user->hin,"\% hin\n");
        fprintf(f,  "%-24.16f   %s",user->hout,"\% hout\n");
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
 
    fclaw2d_domain_barrier (glob->domain);  /* redundant?  */
    RADIALDAM_SETPROB();
}


void radialdam_patch_setup(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *patch,
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

    double *xd,*yd,*zd;
    double *xp,*yp,*zp, *area;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *xnormals,*ynormals,*xtangents,*ytangents;
    double *surfnormals,*edgelengths,*curvature;
    fclaw2d_clawpatch_metric_data2(glob,patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,&edgelengths,
                                   &curvature);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);
    

    const user_options_t* user = radialdam_get_options(glob);
    if (user->claw_version == 4)
        USER46_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                               &dx,&dy,&maux,aux,
                               xnormals,xtangents,
                               ynormals,ytangents,
                               surfnormals,area);

    else if (user->claw_version == 5)
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                              &dx,&dy,&maux,aux,
                              xnormals,xtangents,
                              ynormals,ytangents,
                              surfnormals,area);
    
}