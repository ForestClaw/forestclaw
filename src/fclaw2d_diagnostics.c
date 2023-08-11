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

#include <fclaw2d_diagnostics.h>

#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw2d_options.h>

#include <fclaw_gauges.h>
#include <fclaw_pointer_map.h>

static
fclaw2d_diagnostics_vtable_t* diagnostics_vt_new()
{
    return (fclaw2d_diagnostics_vtable_t*) FCLAW_ALLOC_ZERO (fclaw2d_diagnostics_vtable_t, 1);
}

static
void diagnostics_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

/* ---------------------------------------------------
    Public interface
    ------------------------------------------------ */

fclaw2d_diagnostics_vtable_t* fclaw2d_diagnostics_vt(fclaw_global_t* glob)
{
	fclaw2d_diagnostics_vtable_t* diagnostics_vt = (fclaw2d_diagnostics_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fclaw2d_diagnostics");
	FCLAW_ASSERT(diagnostics_vt != NULL);
	FCLAW_ASSERT(diagnostics_vt->is_set != 0);
    return diagnostics_vt;
}

/* global_maximum is in forestclaw2d.c */
double fclaw2d_domain_global_minimum (fclaw_domain_t* domain, double d)
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
void fclaw2d_diagnostics_vtable_initialize(fclaw_global_t* glob)
{

    fclaw2d_diagnostics_vtable_t *diag_vt = diagnostics_vt_new();

    /* Function pointers all set to zero, since diag_vt has static storage */

    diag_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fclaw2d_diagnostics") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fclaw2d_diagnostics", diag_vt, diagnostics_vt_destroy);
}

static void acc_destroy(void* acc)
{
    FCLAW_FREE(acc);
}

void fclaw2d_diagnostics_initialize(fclaw_global_t *glob)
{
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt(glob);

    fclaw2d_diagnostics_accumulator_t *acc = FCLAW_ALLOC (fclaw2d_diagnostics_accumulator_t, 1);
    fclaw_global_attribute_store(glob, "acc", acc, acc_destroy);
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    /* Return an error accumulator */
    if (diag_vt->patch_init_diagnostics != NULL)
        diag_vt->patch_init_diagnostics(glob,&acc->patch_accumulator);

    /* gauges */
    if (diag_vt->gauges_init_diagnostics != NULL)
        diag_vt->gauges_init_diagnostics(glob,&acc->gauge_accumulator);

    /* rays */
    if (diag_vt->ray_init_diagnostics != NULL)
        diag_vt->ray_init_diagnostics(glob,&acc->ray_accumulator);

    /* solvers */
    if (diag_vt->solver_init_diagnostics != NULL)
        diag_vt->solver_init_diagnostics(glob,&acc->solver_accumulator);

    /* User diagnostics */
    if (fclaw_opt->run_user_diagnostics != 0 && diag_vt->user_init_diagnostics != NULL)
        diag_vt->user_init_diagnostics(glob,&acc->user_accumulator);
}


/* Collect statistics in the accumulator */
void fclaw2d_diagnostics_gather(fclaw_global_t *glob,
                                int init_flag)
{
    fclaw2d_diagnostics_accumulator_t *acc = fclaw_global_get_attribute(glob, "acc");
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt(glob);

    /* -----------------------------------------------------
       Compute diagnostics on all local patches
       ----------------------------------------------------- */
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_DIAGNOSTICS]);


    /* Patch diagnostics */
    if (diag_vt->patch_compute_diagnostics != NULL)
        diag_vt->patch_compute_diagnostics(glob,acc->patch_accumulator);

    /* Gauges diagnostics */
    if (diag_vt->gauges_compute_diagnostics != NULL)
        diag_vt->gauges_compute_diagnostics(glob,acc->gauge_accumulator);

    /* Rays diagnostics */
    if (diag_vt->ray_compute_diagnostics != NULL)
        diag_vt->ray_compute_diagnostics(glob,acc->ray_accumulator);

    /* Solver diagnostics */
    if (diag_vt->solver_compute_diagnostics != NULL)
        diag_vt->solver_compute_diagnostics(glob,acc->solver_accumulator);

    /* User diagnostics */
    if (fclaw_opt->run_user_diagnostics != 0 && diag_vt->user_compute_diagnostics != NULL)
        diag_vt->user_compute_diagnostics(glob,acc->user_accumulator);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_DIAGNOSTICS]);


    /* ---------------------------------------------------------
       Gather all of the statistics (requires communication)
       --------------------------------------------------------- */
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_DIAGNOSTICS_COMM]);

    if (diag_vt->patch_gather_diagnostics != NULL)
        diag_vt->patch_gather_diagnostics(glob,acc->patch_accumulator,init_flag);

    if (diag_vt->gauges_gather_diagnostics != NULL)
        diag_vt->gauges_gather_diagnostics(glob,acc->gauge_accumulator,init_flag);

    if (diag_vt->ray_gather_diagnostics != NULL)
        diag_vt->ray_gather_diagnostics(glob,acc->ray_accumulator,init_flag);

    if (diag_vt->solver_gather_diagnostics != NULL)
        diag_vt->solver_gather_diagnostics(glob,acc->solver_accumulator,init_flag);

    if (fclaw_opt->run_user_diagnostics != 0 && diag_vt->user_gather_diagnostics != NULL)
        diag_vt->user_gather_diagnostics(glob,acc->user_accumulator,init_flag);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_DIAGNOSTICS_COMM]);

    fclaw2d_diagnostics_reset(glob);
}

void fclaw2d_diagnostics_reset(fclaw_global_t *glob)
{
    fclaw2d_diagnostics_accumulator_t *acc = fclaw_global_get_attribute(glob, "acc");
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt(glob);

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_DIAGNOSTICS]);

    if (diag_vt->patch_reset_diagnostics != NULL)
        diag_vt->patch_reset_diagnostics(glob,acc->patch_accumulator);

    if (diag_vt->gauges_reset_diagnostics != NULL)
        diag_vt->gauges_reset_diagnostics(glob,acc->gauge_accumulator);

    if (diag_vt->ray_reset_diagnostics != NULL)
        diag_vt->ray_reset_diagnostics(glob,acc->ray_accumulator);

    if (diag_vt->solver_reset_diagnostics != NULL)
        diag_vt->solver_reset_diagnostics(glob,acc->solver_accumulator);

    if (fclaw_opt->run_user_diagnostics != 0 && diag_vt->user_reset_diagnostics != NULL)
        diag_vt->user_reset_diagnostics(glob,acc->user_accumulator);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_DIAGNOSTICS]);
}

void fclaw2d_diagnostics_finalize(fclaw_global_t *glob)
{
    fclaw2d_diagnostics_accumulator_t *acc = fclaw_global_get_attribute(glob, "acc");
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt(glob);

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_DIAGNOSTICS]);

    if (diag_vt->patch_finalize_diagnostics != NULL)
        diag_vt->patch_finalize_diagnostics(glob,&acc->patch_accumulator);

    if (diag_vt->gauges_finalize_diagnostics != NULL)
        diag_vt->gauges_finalize_diagnostics(glob,&acc->gauge_accumulator);

    if (diag_vt->ray_finalize_diagnostics != NULL)
        diag_vt->ray_finalize_diagnostics(glob,&acc->ray_accumulator);

    if (diag_vt->solver_finalize_diagnostics != NULL)
        diag_vt->solver_finalize_diagnostics(glob,&acc->solver_accumulator);

    if (fclaw_opt->run_user_diagnostics != 0 && diag_vt->user_finalize_diagnostics != NULL)
        diag_vt->user_finalize_diagnostics(glob,&acc->user_accumulator);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_DIAGNOSTICS]);
}

