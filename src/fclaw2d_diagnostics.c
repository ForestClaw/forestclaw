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


/* -----------------------------------------------------------------
   Main routine.

   Note that the check for whether the user has specified diagnostics
   to run is done here, not in fclaw2d_run.cpp
   ---------------------------------------------------------------- */

void* fclaw2d_diagnostics_initialize(fclaw2d_domain_t *domain, int init_flag)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_vtable_t *vt = fclaw2d_vt();

    /* Return an error accumulator */
    return vt.init_diagnostics(domain,init_flag);
}

void fclaw2d_diagnostics_gather(fclaw2d_domain_t *domain,
                                int init_flag,
                                void* gather_accumulator)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_vtable_t *vt = fclaw2d_vt();


    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

    fclaw2d_domain_iterate_patches(domain, cb_fclaw2d_patch_compute_diagnostics,
                                   gather_accumulator);

    vt.gather_diagnostics(domain,gather_accumulator);

    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_DIAGNOSTICS]);

}
