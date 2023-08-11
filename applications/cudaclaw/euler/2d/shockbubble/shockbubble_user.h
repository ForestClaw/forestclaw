/*
  Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef SHOCKBUBBLE_USER_H
#define SHOCKBUBBLE_USER_H

#include <fclaw2d_include_all.h>


#include <fclaw2d_include_all.h>

#include <fc2d_cudaclaw.h>
#include <cudaclaw_user_fort.h>

#include <fc2d_cudaclaw_cuda.h>
#include <fc2d_cudaclaw_options.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <cudaclaw_user_fort.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct user_options
{
    int idisc;
    double gamma;
    double x0;
    double y0;
    double r0;
    double rhoin;
    double pinf;

    int claw_version;
    int cuda;

    int is_registered;
} user_options_t;


void shockbubble_link_solvers(fclaw_global_t *glob);

user_options_t* shockbubble_options_register (fclaw_app_t * app,
                                          const char *configfile);

void shockbubble_options_store (fclaw_global_t* glob, user_options_t* user);

user_options_t* shockbubble_get_options(fclaw_global_t* glob);

/* ------------------------------- Fortran code --------------------------------------- */

#if 0
#define SETPROB FCLAW_F77_FUNC(setprob, SETPROB)
void SETPROB(const double *gamma, const double* x0, const double* y0,
                         const double* r0, const double* rhoin,
                         const double* pinf, const int* idisc);
#endif



/* ---------------------------------- Cuda code --------------------------------------- */
void setprob_cuda();

void shockbubble_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);

void shockbubble_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2);


#ifdef __cplusplus
}
#endif

#endif
