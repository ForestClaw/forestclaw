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

/* -------------------------------------------------------------
   Conservation check
   ------------------------------------------------------------- */

static
void cb_diagnostics_compute_sum(fclaw2d_domain_t *domain,
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
void fclaw2d_check_conservation(fclaw2d_domain_t *domain,
                                int init_flag,
                                fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);

    const amr_options_t *gparms = get_domain_parms(domain);
    int meqn = gparms->meqn;

    double *local_sum = FCLAW_ALLOC_ZERO(double,meqn);
    double *global_sum = FCLAW_ALLOC_ZERO(double,meqn);

    /* This is needed to store the global sum, defined above */
    FCLAW_ASSERT(meqn < FCLAW2D_DIAGNOSTICS_MAX_MEQN);

    /* Accumulate sum for all patches */
    fclaw2d_domain_iterate_patches(domain,cb_diagnostics_compute_sum,(void *) local_sum);


    /* Report results */
    fclaw_global_productionf("Conservation check\n");
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);
    for (int m = 0; m < meqn; m++)
    {
        global_sum[m] = fclaw2d_domain_global_sum (domain, local_sum[m]);
        /* One time setting for initial sum */
        if (init_flag)
        {
            global_sum0[m] = global_sum[m];
        }
        fclaw_global_essentialf("sum[%d] =  %24.16e  %24.16e\n",m,global_sum[m],
                                 fabs(global_sum[m]-global_sum0[m]));
    }
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_start (&ddata->timers[running]);
    }

    fclaw_global_productionf("\n");
    FCLAW_FREE(local_sum);
    FCLAW_FREE(global_sum);

}


/* -------------------------------------------------------------
   Compute total area (needed to normalize errors)

   -- At some point, a option may be provided so that the
      user can supply their own calculation of the total area
   ------------------------------------------------------------- */

static
void cb_diagnostics_patch_area(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               void *user)
{
    double *sum = (double*) user;
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    double *area;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(domain,this_patch);


    *sum += FCLAW2D_FORT_COMPUTE_PATCH_AREA(&mx,&my,&mbc,&dx,&dy,area);
}

static
double fclaw2d_compute_total_area(fclaw2d_domain_t *domain,
                                  fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);

    double local_sum;
    double global_sum;

    /* ----------------------------------------------------
       Accumulate area for all patches on this processor
       ---------------------------------------------------- */
    local_sum = 0;
    fclaw2d_domain_iterate_patches(domain,cb_diagnostics_patch_area,
                                   (void *) &local_sum);


    /* ---------------------------------------
       Compute global sum to get total area
       --------------------------------------- */
    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    /* Do an all gather to get total sum */
    global_sum = fclaw2d_domain_global_sum (domain, local_sum);

    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_start (&ddata->timers[running]);
    }

    return global_sum;
}


/* -------------------------------------------------------------
   Compute errors (1-norm, 2-norm, inf-norm)
   ------------------------------------------------------------- */

typedef struct {
    double* local_error;  /* meqn x 3 array of errors on a patch */
#if 0
    double* error;  /* Grid of errors, computed in user-defined routine */
#endif
} error_info_t;

static
void cb_diagnostics_compute_error(fclaw2d_domain_t *domain,
                                  fclaw2d_patch_t *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx,
                                  void *user)
{
    error_info_t *error_data = (error_info_t*) user;
    double *area, *error;
    double *q;
    int mx, my, mbc, meqn;
    double xlower,ylower,dx,dy;
    fclaw2d_vtable_t vt;

    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(domain,this_patch); /* Might be NULL */
    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    error = fclaw2d_clawpatch_get_error(domain,this_patch);

    vt.compute_patch_error(domain,this_patch,this_block_idx,this_patch_idx,
                           error);

    /* Accumulate sums and maximums needed to compute errors */
    vt.fort_compute_error_norm(&mx, &my, &mbc, &meqn, &dx,&dy, area,
                               error, error_data->local_error);
#if 0
    FCLAW2D_FORT_COMPUTE_ERROR_NORM(&mx, &my, &mbc, &meqn, &dx,&dy, area,
                                    error, error_data->local_error);
#endif
}

static
void fclaw2d_compute_total_error(fclaw2d_domain_t *domain,
                                 fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    const amr_options_t *gparms = get_domain_parms(domain);

    double total_area;
    double *error_norm;
    int meqn;
    error_info_t error_data;

    meqn = gparms->meqn;
    error_norm = FCLAW_ALLOC_ZERO(double,3*meqn);

    /* Allocate memory for 1-norm, 2-norm, and inf-norm errors */
    error_data.local_error  = FCLAW_ALLOC_ZERO(double,3*meqn);

    /* -------------------------------------------------
       Compute error for all patches on this processor
       ------------------------------------------------- */
    fclaw2d_domain_iterate_patches(domain, cb_diagnostics_compute_error,
                                   (void *) &error_data);

    /* -------------------------------------------------
       Accumulate errors on all processors
       ------------------------------------------------- */
    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    for (int m = 0; m < meqn; m++)
    {
        int i1 = m;            /* 1-norm */
        int i2 = meqn + m;     /* 2-norm */
        int iinf = 2*meqn + m; /* inf-norm */
        error_norm[i1]   = fclaw2d_domain_global_sum     (domain, error_data.local_error[i1]);
        error_norm[i2]   = fclaw2d_domain_global_sum     (domain, error_data.local_error[i2]);
        error_norm[iinf] = fclaw2d_domain_global_maximum (domain, error_data.local_error[iinf]);
    }

    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);
    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_start (&ddata->timers[running]);
    }

    /* -----------------------------
       Normalize errors by area
       ----------------------------- */
    total_area = fclaw2d_compute_total_area(domain,
                                            running);
    FCLAW_ASSERT(total_area != 0);

    for(int m = 0; m < meqn; m++)
    {
        int i1 = m;            /* 1-norm */
        int i2 = meqn + m;     /* 2-norm */
        int iinf = 2*meqn + m; /* inf-norm */
        error_norm[i1] = error_norm[i1]/total_area;
        error_norm[i2] = sqrt(error_norm[i2]/total_area);
        fclaw_global_essentialf("error[%d] =  %8.4e  %8.4e %8.4e\n",m,
                                 error_norm[i1], error_norm[i2],error_norm[iinf]);
    }
    FCLAW_FREE(error_norm);
    FCLAW_FREE(error_data.local_error);
}

/* -----------------------------------------------------------------
   Main routine.

   Note that the check for whether the user has specified diagnostics
   to run is done here, not in fclaw2d_run.cpp
   ---------------------------------------------------------------- */

void fclaw2d_diagnostics_run(fclaw2d_domain_t *domain, int init_flag)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_vtable_t vt;

    if (gparms->compute_error)
    {
        fclaw2d_compute_total_error(domain,FCLAW2D_TIMER_DIAGNOSTICS);
    }

    if (gparms->conservation_check)
    {
        fclaw2d_check_conservation(domain,init_flag,FCLAW2D_TIMER_DIAGNOSTICS);
    }

    if (gparms->run_user_diagnostics)
    {
        /* The user could also do this inside of their diagnostic
           function, but it is supplied here to be nice */
        double t = fclaw2d_domain_get_time(domain);

        vt = fclaw2d_get_vtable(domain);

        FCLAW_ASSERT(vt.run_user_diagnostics != NULL);
        vt.run_user_diagnostics(domain,t);
    }

    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);
}
