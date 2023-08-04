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

#include "swirl_user.h"

#include <fclaw3d_metric.h>

void swirl_problem_setup(fclaw2d_global_t *glob)
{
    const user_options_t* user = swirl_get_options(glob);
    fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",fclaw_opt->manifold,"\% manifold\n");

        /* For five-patch mapping */
        fprintf(f,  "%-24.4f   %s",user->alpha,"\% alpha\n");

        /* For Bilinear mapping */
        fprintf(f,  "%-24.4f   %s",user->center[0],"\% center_x\n");
        fprintf(f,  "%-24.4f   %s",user->center[1],"\% center_y\n");
        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);
    SETPROB();
}

static
void swirl_patch_setup(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *patch,
                          int blockno,
                          int patchno)
{
    int mx,my,mz, mbc;
    double xlower,ylower,zlower, dx,dy, dz;
    fclaw3dx_clawpatch_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    double *xd,*yd,*zd,*volume,*faceareas;
    double *xp,*yp,*zp;
    fclaw3d_clawpatch_mesh_data(glob,patch,&xp,&yp,&zp,
                                &xd,&yd,&zd,&volume,&faceareas);

    double *xrot, *yrot, *zrot;
    fclaw3d_metric_patch_basis(glob,patch,&xrot,&yrot,&zrot);

    int maux;
    double *aux;
    fclaw3dx_clawpatch_aux_data(glob,patch,&aux,&maux);

    SWIRL_SETAUX_MANIFOLD(&mbc,&mx,&my,&mz, &xlower,&ylower,&zlower,
                          &dx,&dy,&dz,&maux,aux,&blockno,
                          xrot,yrot,zrot,volume,faceareas);
}


void swirl_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt(glob);
    fclaw_vt->problem_setup = swirl_problem_setup;

    /* example of how to set up a user defined criteria */
    fclaw3dx_clawpatch_vtable_t *clawpatch_vt = fclaw3dx_clawpatch_vt(glob);
    clawpatch_vt->d3->fort_user_exceeds_threshold = &FCLAW3DX_USER_EXCEEDS_TH;

    fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

    const user_options_t* user = swirl_get_options(glob);
    if (user->claw_version == 4)
    {
        fc3d_clawpack46_vtable_t *clawpack46_vt = fc3d_clawpack46_vt(glob);        

        if (fclaw_opt->manifold)
        {
            /* This calls a manifold version of setaux */
            fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
            patch_vt->setup = swirl_patch_setup;
        }
        else
        {
            /* This calls a setaux routine */
            clawpack46_vt->fort_setaux     = &CLAWPACK46_SETAUX;  
        }

        clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->fort_rpn3      = &CLAWPACK46_RPN3;
        clawpack46_vt->fort_rpt3      = &CLAWPACK46_RPT3;
        clawpack46_vt->fort_rptt3      = &CLAWPACK46_RPTT3;
    }
    else if (user->claw_version == 5)
    {
        printf("swirl_user.cpp : Example not implemented for Claw version 5.\n");
        exit(0);
    }
}





