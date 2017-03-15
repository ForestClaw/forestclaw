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

#ifndef FCLAW2D_DIAGNOSTICS_H
#define FCLAW2D_DIAGNOSTICS_H

#include "forestclaw2d.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fclaw2d_diagnostics_vtable  fclaw2d_diagnostics_vtable_t;

typedef struct
{
    void* patch_accumulator;
    void* solver_accumulator;
    void* user_accumulator;
} fclaw2d_diagnostics_accumulator_t;

/* Diagnostic information */
typedef void (*fclaw2d_diagnostics_initialize_t)(fclaw2d_global_t *glob,
                                                 void** acc);

typedef void (*fclaw2d_diagnostics_compute_t)(fclaw2d_global_t *glob,
                                              void* acc);

typedef void (*fclaw2d_diagnostics_gather_t)(fclaw2d_global_t *glob,
                                             void* acc,
                                             int init_flag);

typedef void (*fclaw2d_diagnostics_reset_t)(fclaw2d_global_t *glob,
                                            void* acc);

typedef void (*fclaw2d_diagnostics_finalize_t)(fclaw2d_global_t *glob,
                                               void** acc);

fclaw2d_diagnostics_vtable_t* fclaw2d_diagnostics_vt();


struct fclaw2d_diagnostics_vtable
{
    /* patch diagnostic functions (error, conservation, area, etc) */
    fclaw2d_diagnostics_initialize_t     patch_init_diagnostics;
    fclaw2d_diagnostics_compute_t        patch_compute_diagnostics;
    fclaw2d_diagnostics_gather_t         patch_gather_diagnostics;
    fclaw2d_diagnostics_reset_t          patch_reset_diagnostics;
    fclaw2d_diagnostics_finalize_t       patch_finalize_diagnostics;

    /* solver diagnostic functions (gauges, fgmax, and so on) */
    fclaw2d_diagnostics_initialize_t     solver_init_diagnostics;
    fclaw2d_diagnostics_compute_t        solver_compute_diagnostics;
    fclaw2d_diagnostics_gather_t         solver_gather_diagnostics;
    fclaw2d_diagnostics_reset_t          solver_reset_diagnostics;
    fclaw2d_diagnostics_finalize_t       solver_finalize_diagnostics;

    /* user defined diagnostics */
    fclaw2d_diagnostics_initialize_t     user_init_diagnostics;
    fclaw2d_diagnostics_compute_t        user_compute_diagnostics;
    fclaw2d_diagnostics_gather_t         user_gather_diagnostics;
    fclaw2d_diagnostics_reset_t          user_reset_diagnostics;
    fclaw2d_diagnostics_finalize_t       user_finalize_diagnostics;

};

void fclaw2d_diagnostics_vtable_init();

/* See forestclaw2d.h for the maximum version of this function */
double fclaw2d_domain_global_minimum (fclaw2d_domain_t* domain, double d);

void fclaw2d_diagnostics_initialize(fclaw2d_global_t *glob,
                                    fclaw2d_diagnostics_accumulator_t* acc);

void fclaw2d_diagnostics_gather(fclaw2d_global_t *glob,
                                fclaw2d_diagnostics_accumulator_t* acc,
                                int init_flag);

void fclaw2d_diagnostics_reset(fclaw2d_global_t *glob,
                               fclaw2d_diagnostics_accumulator_t* acc);

void fclaw2d_diagnostics_finalize(fclaw2d_global_t *glob,
                                  fclaw2d_diagnostics_accumulator_t* acc);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
