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

#include "transport_user.h"

void transport_b4step2_manifold(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *patch,
                                int blockno,
                                int patchno,
                                double t, double dt,
                                int claw_version)

{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw2d_clawpatch_metric_vector(glob,patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    double *xp, *yp, *zp, *xd, *yd, *zd, *area;
    fclaw2d_clawpatch_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);

    // const user_options_t* user_opt = sphere_get_options(glob);
    if (claw_version == 4)
    {
        TRANSPORT46_SET_VELOCITY(&blockno, &mx, &my, &mbc,
                                 &dx, &dy, &xlower, &ylower,
                                 &t, xp, yp, zp, 
                                 xnormals,ynormals, surfnormals,
                                 aux,&maux);
    }
    else
    {
        TRANSPORT5_SET_VELOCITY(&blockno, &mx, &my, &mbc,
                                &dx, &dy, &xlower, &ylower,
                                &t, xp, yp, zp, 
                                xnormals,ynormals, surfnormals,
                                aux,&maux);
    }    
}