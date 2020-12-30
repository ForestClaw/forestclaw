/*
Copyright (c) 2012-2020 Carsten Burstedde, Donna Calhoun
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

#include "heat_diagnostics.h"

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>


#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_diagnostics.h>

void heat_diagnostics_initialize(fclaw2d_global_t *glob,
                                   void **acc_patch)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
              fclaw2d_clawpatch_get_options(glob);

    heat_error_info_t *error_data;

    int mfields = clawpatch_opt->rhs_fields;

    error_data = FCLAW_ALLOC(heat_error_info_t,1);

    /* Allocate memory for 1-norm, 2-norm, and inf-norm errors */
    error_data->local_error  = FCLAW_ALLOC_ZERO(double,3*mfields);
    error_data->global_error  = FCLAW_ALLOC_ZERO(double,3*mfields);
    error_data->mass   = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->mass0  = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->rhs   = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->boundary   = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->area = 0;
    error_data->c_kahan = FCLAW_ALLOC_ZERO(double,mfields);   

    *acc_patch = error_data;

}

void heat_diagnostics_reset(fclaw2d_global_t *glob,
                              void* patch_acc)
{
    heat_error_info_t *error_data = (heat_error_info_t*) patch_acc;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    int mfields = clawpatch_opt->rhs_fields;

    for(int m = 0; m < mfields; m++)
    {
        int i1 = m;            /* 1-norm */
        int i2 = mfields + m;     /* 2-norm */
        int iinf = 2*mfields + m; /* inf-norm */
        error_data->local_error[i1] = 0;
        error_data->local_error[i2] = 0;
        error_data->local_error[iinf] = 0;
        error_data->mass[m]= 0;        

        error_data->rhs[m] = 0;
        error_data->boundary[m] = 0;
        error_data->c_kahan[m] = 0;
    }
    error_data->area = 0;
}

static
void heat_compute(fclaw2d_domain_t *domain,
                    fclaw2d_patch_t *patch,
                    int blockno,
                    int patchno,
                    void* user) //void *patch_acc)
{
    fclaw2d_global_iterate_t *s = (fclaw2d_global_iterate_t *) user;
    heat_error_info_t *error_data = (heat_error_info_t*) s->user; 

    /* Accumulate area for final computation of error */
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(s->glob,patch,&mx,&my,&mbc,&xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    double *area = fclaw2d_clawpatch_get_area(s->glob,patch);  
    FCLAW_ASSERT(clawpatch_vt->fort_compute_patch_area != NULL);
    error_data->area += clawpatch_vt->fort_compute_patch_area(&mx,&my,&mbc,&dx,&dy,area);

    /* Compute error */
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(s->glob);
    if (fclaw_opt->compute_error)
    {
        clawpatch_vt->compute_error(s->glob, patch, blockno, patchno, error_data);
    }

    if (fclaw_opt->conservation_check)
    {
        clawpatch_vt->conservation_check(s->glob, patch, blockno, patchno, error_data);
    }
}

void heat_diagnostics_compute(fclaw2d_global_t* glob,
                                           void* patch_acc)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int check = fclaw_opt->compute_error || fclaw_opt->conservation_check;
    if (!check) return;

    fclaw2d_global_iterate_patches(glob, heat_compute, patch_acc);
}



/* Accumulate the errors computed above */
void heat_diagnostics_gather(fclaw2d_global_t *glob,
                               void* patch_acc,
                               int init_flag)
{
    fclaw2d_domain_t *domain = glob->domain;
    
    heat_error_info_t *error_data = (heat_error_info_t*) patch_acc;
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    
    int mfields = clawpatch_opt->rhs_fields;  /* clawpatch->meqn */

    if (fclaw_opt->compute_error != 0)
    {
        double total_area = fclaw2d_domain_global_sum(domain, error_data->area);
        FCLAW_ASSERT(total_area != 0);

        double *error_norm = FCLAW_ALLOC_ZERO(double,3*mfields);
        for (int m = 0; m < mfields; m++)
        {
            int i1 = m;            /* 1-norm */
            int i2 = mfields + m;     /* 2-norm */
            int i3 = 2*mfields + m; /* inf-norm */

            error_norm[i1]  = fclaw2d_domain_global_sum(domain, error_data->local_error[i1]);
            error_norm[i1] /= total_area;

            error_norm[i2]  = fclaw2d_domain_global_sum(domain, error_data->local_error[i2]);
            error_norm[i2] /= total_area;
            error_norm[i2] = sqrt(error_norm[i2]);

            error_norm[i3] = fclaw2d_domain_global_maximum(domain, 
                                                           error_data->local_error[i3]);

            error_data->global_error[i1] = error_norm[i1];
            error_data->global_error[i2] = error_norm[i2];
            error_data->global_error[i3] = error_norm[i3];

            fclaw_global_essentialf("error[%d] = %16.6e %16.6e %16.6e\n",m,
                                    error_norm[i1], error_norm[i2],error_norm[i3]);
        }
        FCLAW_FREE(error_norm);
    }


    if (fclaw_opt->conservation_check != 0)
    {
        double *total_mass = FCLAW_ALLOC_ZERO(double,mfields);
        for(int m = 0; m < mfields; m++)
        {
            /* Store mass for future checks */
            if (init_flag)
            {
                total_mass[m] = fclaw2d_domain_global_sum(domain, error_data->rhs[m]);
                error_data->mass0[m] = total_mass[m];                
            }
            else
            {
                total_mass[m] = fclaw2d_domain_global_sum(domain, error_data->boundary[m]);
            }
            fclaw_global_essentialf("sum[%d] =  %24.16e  %24.16e\n",m,total_mass[m],
                                    fabs(total_mass[m]-error_data->mass0[m]));
        }
        FCLAW_FREE(total_mass);        
    }
}


void heat_diagnostics_finalize(fclaw2d_global_t *glob,
                                 void** patch_acc)
{
    heat_error_info_t *error_data = *((heat_error_info_t**) patch_acc);
    FCLAW_FREE(error_data->mass);
    FCLAW_FREE(error_data->mass0);
    FCLAW_FREE(error_data->rhs);
    FCLAW_FREE(error_data->boundary);
    FCLAW_FREE(error_data->c_kahan);
    FCLAW_FREE(error_data->local_error);
    FCLAW_FREE(error_data->global_error);
    FCLAW_FREE(error_data);
    *patch_acc = NULL;
}
