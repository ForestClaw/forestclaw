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

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_domain.h>
#include <fclaw_diagnostics.h>

void fclaw_clawpatch_diagnostics_cons_default(fclaw_global_t *glob,
                                              fclaw_patch_t *patch,
                                              int blockno,
                                              int patchno,
                                              void *user)
{
    error_info_t* error_data = (error_info_t*) user;
    double* area = fclaw_clawpatch_get_2d_area(glob,patch);  /* Might be null */

    int meqn;
    double *q; 
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    if(clawpatch_vt->patch_dim == 2)
    {
        FCLAW_ASSERT(clawpatch_vt->d2->fort_conservation_check != NULL);
        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);
        clawpatch_vt->d2->fort_conservation_check(&mx, &my, &mbc, &meqn, &dx,&dy,
                                                  area, q, error_data->mass,
                                                  error_data->c_kahan);
    }
    else
    {
        FCLAW_ASSERT(clawpatch_vt->d3->fort_conservation_check != NULL);
        int mz;
        double zlower, dz;
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                    &xlower,&ylower,&zlower, &dx,&dy, &dz);    
        clawpatch_vt->d3->fort_conservation_check(&mx, &my, &mz, &mbc, &meqn, 
                                                  &dx,&dy,&dz,
                                                  area, q, error_data->mass,
                                                  error_data->c_kahan);
    }
}

void fclaw_clawpatch_diagnostics_error_default(fclaw_global_t *glob,
                                               fclaw_patch_t *patch,
                                               int blockno,
                                               int patchno,
                                               void *user)
{
    error_info_t* error_data = (error_info_t*) user;
    //const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    double *area = fclaw_clawpatch_get_2d_area(glob,patch);  /* Might be null */
    double *q;
    int meqn;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    if (clawpatch_vt->patch_dim == 2 && clawpatch_vt->d2->fort_compute_patch_error != NULL)
    {
        double t = glob->curr_time;
        double* error = fclaw_clawpatch_get_error(glob,patch);
        double* soln = fclaw_clawpatch_get_exactsoln(glob,patch);

        int mx, my, mbc;
        double xlower,ylower,dx,dy;

        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,&xlower,&ylower,&dx,&dy);

        clawpatch_vt->d2->fort_compute_patch_error(&blockno, &mx,&my,&mbc,&meqn,&dx,&dy,
                                              &xlower,&ylower, &t, q, error, soln);

        /* Accumulate sums and maximums needed to compute error norms */
        FCLAW_ASSERT(clawpatch_vt->d2->fort_compute_error_norm != NULL);
        clawpatch_vt->d2->fort_compute_error_norm(&blockno, &mx, &my, &mbc, &meqn, 
                                                  &dx,&dy, area, error,
                                                  error_data->local_error);
    }
    else if (clawpatch_vt->patch_dim == 3 && clawpatch_vt->d3->fort_compute_patch_error != NULL)
    {
        double t = glob->curr_time;
        double* error = fclaw_clawpatch_get_error(glob,patch);
        double* soln = fclaw_clawpatch_get_exactsoln(glob,patch);

        int mx, my, mz, mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, 
                                    &mbc,&xlower,&ylower,&zlower, 
                                    &dx,&dy,&dz);
        clawpatch_vt->d3->fort_compute_patch_error(&blockno, &mx,&my,&mz,
                                               &mbc,&meqn,
                                               &dx,&dy,&dz,
                                               &xlower,&ylower, &zlower,
                                               &t, q, error, soln);

        /* Accumulate sums and maximums needed to compute error norms */
        FCLAW_ASSERT(clawpatch_vt->d3->fort_compute_error_norm != NULL);
        clawpatch_vt->d3->fort_compute_error_norm(&blockno, &mx, &my, &mz,
                                                  &mbc, &meqn, 
                                                  &dx,&dy, &dz, area, error,
                                                  error_data->local_error);


    }
}


