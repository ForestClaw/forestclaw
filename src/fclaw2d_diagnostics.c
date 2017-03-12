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
#include <fclaw2d_vtable.h>


/* global_maximum is in forestclaw2d.c */
double fclaw2d_domain_global_minimum (fclaw2d_domain_t* domain, double d)
{
    double neg_d;
    double maxvalue;
    neg_d = -d;
    maxvalue = fclaw2d_domain_global_maximum(domain,neg_d);
    return -maxvalue;
}


/* -----------------------------------------------------------------
   Main routine.

   Note that the check for whether the user has specified diagnostics
   to run is done here, not in fclaw2d_run.cpp
   ---------------------------------------------------------------- */

void fclaw2d_diagnostics_initialize(fclaw2d_domain_t *domain,
                                    fclaw2d_diagnostics_accumulator_t* acc)
{
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    const amr_options_t *gparms = get_domain_parms(domain);

    /* Return an error accumulator */
    if (fclaw_vt->patch_init_diagnostics != NULL)
    {
        fclaw_vt->patch_init_diagnostics(domain,&acc->patch_accumulator);
    }

    if (fclaw_vt->solver_init_diagnostics != NULL)
    {
        /* Gauges, fgmax, etc */
        fclaw_vt->solver_init_diagnostics(domain,&acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && fclaw_vt->user_init_diagnostics != NULL)
    {
        fclaw_vt->user_init_diagnostics(domain,&acc->user_accumulator);
    }
}


/* Collect statistics in the accumulator */
void fclaw2d_diagnostics_gather(fclaw2d_domain_t *domain,
                                fclaw2d_diagnostics_accumulator_t* acc,
                                int init_flag)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();


    /* -----------------------------------------------------
       Compute diagnostics on all local patches
       ----------------------------------------------------- */
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);


    if (fclaw_vt->patch_compute_diagnostics != NULL)
    {
        fclaw_vt->patch_compute_diagnostics(domain,acc->patch_accumulator);
    }

    if (fclaw_vt->solver_compute_diagnostics != NULL)
    {
        fclaw_vt->solver_compute_diagnostics(domain,acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && fclaw_vt->user_compute_diagnostics != NULL)
    {
        fclaw_vt->user_compute_diagnostics(domain,acc->user_accumulator);
    }

    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);


    /* ---------------------------------------------------------
       Gather all of the statistics (requires communication)
       --------------------------------------------------------- */
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    if (fclaw_vt->patch_gather_diagnostics != NULL)
    {
        fclaw_vt->patch_gather_diagnostics(domain,acc->patch_accumulator,init_flag);
    }

    if (fclaw_vt->solver_gather_diagnostics != NULL)
    {
        fclaw_vt->solver_gather_diagnostics(domain,acc->solver_accumulator,init_flag);
    }

    if (gparms->run_user_diagnostics != 0 && fclaw_vt->user_gather_diagnostics != NULL)
    {
        fclaw_vt->user_gather_diagnostics(domain,acc->user_accumulator,init_flag);
    }

    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    fclaw2d_diagnostics_reset(domain,acc);
}

void fclaw2d_diagnostics_reset(fclaw2d_domain_t *domain,
                               fclaw2d_diagnostics_accumulator_t* acc)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();

    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

    if (fclaw_vt->patch_reset_diagnostics != NULL)
    {
        fclaw_vt->patch_reset_diagnostics(domain,acc->patch_accumulator);
    }

    if (fclaw_vt->solver_reset_diagnostics != NULL)
    {
        fclaw_vt->solver_reset_diagnostics(domain,acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && fclaw_vt->user_reset_diagnostics != NULL)
    {
        fclaw_vt->user_reset_diagnostics(domain,acc->user_accumulator);
    }
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);
}

void fclaw2d_diagnostics_finalize(fclaw2d_domain_t *domain,
                                  fclaw2d_diagnostics_accumulator_t* acc)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();

    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

    if (fclaw_vt->patch_finalize_diagnostics != NULL)
    {
        fclaw_vt->patch_finalize_diagnostics(domain,&acc->patch_accumulator);
    }

    if (fclaw_vt->solver_finalize_diagnostics != NULL)
    {
        fclaw_vt->solver_finalize_diagnostics(domain,&acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && fclaw_vt->user_finalize_diagnostics != NULL)
    {
        fclaw_vt->user_finalize_diagnostics(domain,&acc->user_accumulator);
    }
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);
}
