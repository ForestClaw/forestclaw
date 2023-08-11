/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#include "advection_user.h"

void claw3_advection_patch_setup_manifold(fclaw2d_global_t *glob,
                                          fclaw_patch_t *patch,
                                          int blockno,
                                          int patchno,
                                          int claw_version)
{
    int mx,my,mz, mbc;
    double xlower,ylower,zlower, dx,dy, dz;
    fclaw3d_clawpatch_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    double *xd,*yd,*zd,*volume,*faceareas;
    double *xp,*yp,*zp;
    fclaw3d_clawpatch_mesh_data(glob,patch,&xp,&yp,&zp,
                                &xd,&yd,&zd,&volume,&faceareas);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    if (claw_version == 4)
        CLAW3_SETAUX_MANIFOLD(&mbc,&mx,&my,&mz, &xlower,&ylower,&zlower,
                               &dx,&dy,&dz,&maux,aux,&blockno,
                               xd,yd,zd,xp,yp,zp,volume,faceareas);

    else if (claw_version == 5)
        fclaw_global_essentialf("claw3_patch_setup_manifold : claw 5 not implemented.\n");
#if 0
        USER5_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                              &dx,&dy,&maux,aux,&blockno,
                              xd,yd,zd,area);
#endif
}





