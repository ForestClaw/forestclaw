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
#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 2
#endif

#if REFINE_DIM == 2 && PATCH_DIM == 2

#include <fclaw2d_clawpatch_diagnostics.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#elif REFINE_DIM == 2 && PATCH_DIM == 3

#include <fclaw3dx_clawpatch_diagnostics.h>

#include <fclaw3dx_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <_fclaw2d_to_fclaw3dx.h>
#include <_fclaw2d_to_fclaw3d.h>

#endif

#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_diagnostics.h>

void fclaw2d_clawpatch_diagnostics_cons_default(fclaw2d_global_t *glob,
                                                fclaw2d_patch_t *patch,
                                                int blockno,
                                                int patchno,
                                                void *user)
{
    error_info_t* error_data = (error_info_t*) user;
    double* area = fclaw2d_clawpatch_get_area(glob,patch);  /* Might be null */

    int meqn;
    double *q; 
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    int mx, my, mbc;
    double xlower,ylower,dx,dy;
#if PATCH_DIM == 2
    FCLAW_ASSERT(clawpatch_vt->d2->fort_conservation_check != NULL);
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    clawpatch_vt->d2->fort_conservation_check(&mx, &my, &mbc, &meqn, &dx,&dy,
                                              area, q, error_data->mass,
                                              error_data->c_kahan);
#else
    FCLAW_ASSERT(clawpatch_vt->d3->fort_conservation_check != NULL);
    int mz;
    double zlower, dz;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);    
    clawpatch_vt->d3->fort_conservation_check(&mx, &my, &mz, &mbc, &meqn, 
                                              &dx,&dy,&dz,
                                              area, q, error_data->mass,
                                              error_data->c_kahan);

#endif

}

void fclaw2d_clawpatch_diagnostics_error_default(fclaw2d_global_t *glob,
                                                 fclaw2d_patch_t *patch,
                                                 int blockno,
                                                 int patchno,
                                                 void *user)
{
    error_info_t* error_data = (error_info_t*) user;
    //const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    double *area = fclaw2d_clawpatch_get_area(glob,patch);  /* Might be null */
    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

#if PATCH_DIM == 2        
    if (clawpatch_vt->d2->fort_compute_patch_error != NULL)
#else
    if (clawpatch_vt->d3->fort_compute_patch_error != NULL)
#endif
    {
        double t = glob->curr_time;
        double* error = fclaw2d_clawpatch_get_error(glob,patch);
        double* soln = fclaw2d_clawpatch_get_exactsoln(glob,patch);

        int mx, my, mbc;
        double xlower,ylower,dx,dy;
#if PATCH_DIM == 2        
        fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,&xlower,&ylower,&dx,&dy);

        clawpatch_vt->d2->fort_compute_patch_error(&blockno, &mx,&my,&mbc,&meqn,&dx,&dy,
                                              &xlower,&ylower, &t, q, error, soln);

        /* Accumulate sums and maximums needed to compute error norms */
        FCLAW_ASSERT(clawpatch_vt->d2->fort_compute_error_norm != NULL);
        clawpatch_vt->d2->fort_compute_error_norm(&blockno, &mx, &my, &mbc, &meqn, 
                                                  &dx,&dy, area, error,
                                                  error_data->local_error);
#elif PATCH_DIM == 3
        int mz;
        double zlower, dz;
        fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mz, 
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


#endif
    }
}


