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

#include <fclaw2d_domain.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_diagnostics_fort.h>
#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_vtable.h>


/* Avoid problem of deleting memory at the end of the run */
#define FCLAW2D_DIAGNOSTICS_MAX_MEQN 20
static double global_sum0[FCLAW2D_DIAGNOSTICS_MAX_MEQN];


/* global_maximum is in forestclaw2d.c */
double fclaw2d_domain_global_minimum (fclaw2d_domain_t* domain, double d)
{
    double neg_d;
    double maxvalue;
    neg_d = -d;
    maxvalue = fclaw2d_domain_global_maximum(domain,neg_d);
    return -maxvalue;
}

static
void fclaw2d_diagnostics_compute_sum(fclaw2d_domain_t *domain,
                                     fclaw2d_patch_t *this_patch,
                                     int this_block_idx,
                                     int this_patch_idx,
                                     void *user)
{
    double *sum = (double*) user;
    double *area;
    double *q;
    int mx, my, mbc, meqn;
    double xlower,ylower,dx,dy;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(domain,this_patch);
    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    FCLAW2D_FORT_CONSERVATION_CHECK(&mx, &my, &mbc, &meqn, &dx,&dy, area, q,sum);
}

static
void fclaw2d_check_conservation(fclaw2d_domain_t *domain, const double t, int init_flag)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int meqn = gparms->meqn;

    double *local_sum = FCLAW_ALLOC_ZERO(double,meqn);
    double *global_sum = FCLAW_ALLOC_ZERO(double,meqn);

    FCLAW_ASSERT(meqn < FCLAW2D_DIAGNOSTICS_MAX_MEQN);

    /* Accumulate sum for all patches */
    fclaw2d_domain_iterate_patches(domain,fclaw2d_diagnostics_compute_sum,(void *) local_sum);


    /* Report results */
    fclaw_global_productionf("Conservation check\n");
    for (int i = 0; i < meqn; i++)
    {
        global_sum[i] = fclaw2d_domain_global_sum (domain, local_sum[i]);
        /* One time setting for initial sum */
        if (init_flag)
        {
            global_sum0[i] = global_sum[i];
        }
        fclaw_global_productionf("sum[%d] =  %24.16e  %24.16e\n",i,global_sum[i],
                                 fabs(global_sum[i]-global_sum0[i]));
    }

    fclaw_global_productionf("\n");
    FCLAW_FREE(local_sum);
    FCLAW_FREE(global_sum);

}

/* -----------------------------------------------------------------
   Main routine.

   Note that the check for whether the user has specified diagnostics
   to run is done here, not in fclaw2d_run.cpp
   ---------------------------------------------------------------- */

void fclaw2d_run_diagnostics(fclaw2d_domain_t *domain, int init_flag)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_vtable_t vt;
    double t;

    t = fclaw2d_domain_get_time(domain);

    if (gparms->run_user_diagnostics)
    {
        vt = fclaw2d_get_vtable(domain);

        FCLAW_ASSERT(vt.run_diagnostics != NULL);
        vt.run_diagnostics(domain,t);
    }

    if (gparms->conservation_check)
    {
        fclaw2d_check_conservation(domain,t,init_flag);
    }

}
