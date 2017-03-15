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

static fclaw2d_diagnostics_vtable_t s_diag_vt;

static
fclaw2d_diagnostics_vtable_t* diagnostics_vt()
{
    return &s_diag_vt;
}

fclaw2d_diagnostics_vtable_t* fclaw2d_diagnostics_vt()
{
    return diagnostics_vt();
}



/* global_maximum is in forestclaw2d.c */
double fclaw2d_domain_global_minimum (fclaw2d_domain_t* domain, double d)
{
    double neg_d;
    double maxvalue;
    neg_d = -d;
    maxvalue = fclaw2d_domain_global_maximum(domain,neg_d);
    return -maxvalue;
}

void fclaw2d_diagnostics_vtable_init()
{

    fclaw2d_diagnostics_vtable_t *diag_vt = diagnostics_vt();

    /* Patch diagnostics (conservation; errors) */
    diag_vt->patch_init_diagnostics       = NULL;  /* allocate memory; set accumulators to 0 */
    diag_vt->patch_compute_diagnostics    = NULL;  /* Update local accumulators */
    diag_vt->patch_gather_diagnostics     = NULL;  /* Gather accumulators from procs; report */
    diag_vt->patch_reset_diagnostics      = NULL;  /* Reset accumulators (e.g. to 0) */
    diag_vt->patch_finalize_diagnostics   = NULL;  /* De-allocate memory */

    /* Solver diagnostics (gauges, fixed-grids (fgmax) */
    diag_vt->solver_init_diagnostics      = NULL;
    diag_vt->solver_compute_diagnostics    = NULL;
    diag_vt->solver_gather_diagnostics    = NULL;
    diag_vt->solver_reset_diagnostics     = NULL;
    diag_vt->solver_finalize_diagnostics  = NULL;

    /* User defined diagnostics (either patch based or otherwise) */
    diag_vt->user_init_diagnostics       = NULL;
    diag_vt->user_compute_diagnostics    = NULL;
    diag_vt->user_gather_diagnostics     = NULL;
    diag_vt->user_reset_diagnostics      = NULL;
    diag_vt->user_finalize_diagnostics   = NULL;
}


/* -----------------------------------------------------------------
   Main routine.

   Note that the check for whether the user has specified diagnostics
   to run is done here, not in fclaw2d_run.cpp
   ---------------------------------------------------------------- */
static
int run_diagnostics(fclaw2d_global_t* glob)
{
    /* Check to see if we should be running any diagnostics */
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);

    int run_diag = gparms->conservation_check ||
                   gparms->run_user_diagnostics ||
                   gparms->compute_error;
    return run_diag;

}



void fclaw2d_diagnostics_initialize(fclaw2d_global_t *glob)
{
    fclaw2d_diagnostics_vtable_t *diag_vt = diagnostics_vt();
    fclaw2d_diagnostics_accumulator_t *acc = glob->acc;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
#if 0
    int run_diag = run_diagnostics(glob);
    if (!run_diag)
    {
        return;
    }
#endif
    /* Return an error accumulator */
    if (diag_vt->patch_init_diagnostics != NULL)
    {
        diag_vt->patch_init_diagnostics(glob,&acc->patch_accumulator);
    }

    if (diag_vt->solver_init_diagnostics != NULL)
    {
        /* Gauges, fgmax, etc */
        diag_vt->solver_init_diagnostics(glob,&acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && diag_vt->user_init_diagnostics != NULL)
    {
        diag_vt->user_init_diagnostics(glob,&acc->user_accumulator);
    }
}


/* Collect statistics in the accumulator */
void fclaw2d_diagnostics_gather(fclaw2d_global_t *glob,
                                int init_flag)
{
    fclaw2d_diagnostics_accumulator_t *acc = glob->acc;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    fclaw2d_diagnostics_vtable_t *diag_vt = diagnostics_vt();

#if 0
    int run_diag = run_diagnostics(glob);
    if (!run_diag)
    {
        return;
    }
#endif

    /* -----------------------------------------------------
       Compute diagnostics on all local patches
       ----------------------------------------------------- */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS]);


    if (diag_vt->patch_compute_diagnostics != NULL)
    {
        diag_vt->patch_compute_diagnostics(glob,acc->patch_accumulator);
    }

    if (diag_vt->solver_compute_diagnostics != NULL)
    {
        diag_vt->solver_compute_diagnostics(glob,acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && diag_vt->user_compute_diagnostics != NULL)
    {
        diag_vt->user_compute_diagnostics(glob,acc->user_accumulator);
    }

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS]);


    /* ---------------------------------------------------------
       Gather all of the statistics (requires communication)
       --------------------------------------------------------- */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    if (diag_vt->patch_gather_diagnostics != NULL)
    {
        diag_vt->patch_gather_diagnostics(glob,acc->patch_accumulator,init_flag);
    }

    if (diag_vt->solver_gather_diagnostics != NULL)
    {
        diag_vt->solver_gather_diagnostics(glob,acc->solver_accumulator,init_flag);
    }

    if (gparms->run_user_diagnostics != 0 && diag_vt->user_gather_diagnostics != NULL)
    {
        diag_vt->user_gather_diagnostics(glob,acc->user_accumulator,init_flag);
    }

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM]);

    fclaw2d_diagnostics_reset(glob);
}

void fclaw2d_diagnostics_reset(fclaw2d_global_t *glob)
{
    fclaw2d_diagnostics_accumulator_t *acc = glob->acc;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    fclaw2d_diagnostics_vtable_t *diag_vt = diagnostics_vt();

    int run_diag = run_diagnostics(glob);
    if (!run_diag)
    {
        return;
    }


    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

    if (diag_vt->patch_reset_diagnostics != NULL)
    {
        diag_vt->patch_reset_diagnostics(glob,acc->patch_accumulator);
    }

    if (diag_vt->solver_reset_diagnostics != NULL)
    {
        diag_vt->solver_reset_diagnostics(glob,acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && diag_vt->user_reset_diagnostics != NULL)
    {
        diag_vt->user_reset_diagnostics(glob,acc->user_accumulator);
    }
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS]);
}

void fclaw2d_diagnostics_finalize(fclaw2d_global_t *glob)
{
    fclaw2d_diagnostics_accumulator_t *acc = glob->acc;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    fclaw2d_diagnostics_vtable_t *diag_vt = diagnostics_vt();

    int run_diag = run_diagnostics(glob);
    if (!run_diag)
    {
        return;
    }


    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

    if (diag_vt->patch_finalize_diagnostics != NULL)
    {
        diag_vt->patch_finalize_diagnostics(glob,&acc->patch_accumulator);
    }

    if (diag_vt->solver_finalize_diagnostics != NULL)
    {
        diag_vt->solver_finalize_diagnostics(glob,&acc->solver_accumulator);
    }

    if (gparms->run_user_diagnostics != 0 && diag_vt->user_finalize_diagnostics != NULL)
    {
        diag_vt->user_finalize_diagnostics(glob,&acc->user_accumulator);
    }
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_DIAGNOSTICS]);
}

