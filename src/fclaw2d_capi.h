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

#ifndef FCLAW2D_CAPI_H
#define FCLAW2D_CAPI_H

/* amr_options.h pulls in sc_{options,obstack,containers}.h and sc.h. */
#if 0
#include "amr_options.h"
#endif
#include "fclaw_options.h"
#include "forestclaw2d.h"
#include "fclaw_base.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* Use as an alternate to GNU feenableexcept */
#ifndef FCLAW_HAVE_FEENABLEEXCEPT
#include <fp_exception_glibc_extension.h>
#endif

#include <fenv.h>
#include <signal.h>

#ifdef FCLAW_HAVE_UNISTD_H
#include <unistd.h>    /* To get process ids */
#endif

/* -----------------------------------------------------------------
  Debug routines
   ---------------------------------------------------------------- */

const amr_options_t* get_domain_parms(fclaw2d_domain_t *domain);
const amr_options_t* fclaw2d_forestclaw_get_options(fclaw2d_domain_t *domain);
void* fclaw2d_domain_get_user_options(fclaw2d_domain_t* domain);
void fclaw2d_domain_set_app(fclaw2d_domain_t* domain,fclaw_app_t* app);
fclaw_app_t* fclaw2d_domain_get_app(fclaw2d_domain_t* domain);

void set_domain_time(fclaw2d_domain_t *domain, double time);
double get_domain_time(fclaw2d_domain_t *domain);

/* Misc. routines */
int num_patches(fclaw2d_domain_t *domain, int level,int include_shadow);
int pow_int(int a, int n);

/* These two are defined in amr_utils.cpp */
void fclaw2d_mpi_debug();

void amrinit(fclaw2d_domain_t **domain);
void amrrun(fclaw2d_domain_t **domain);
void amrreset(fclaw2d_domain_t **domain);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
