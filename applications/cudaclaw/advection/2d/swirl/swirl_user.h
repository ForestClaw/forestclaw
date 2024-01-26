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

#ifndef SWIRL_USER_H
#define SWIRL_USER_H

#include <fclaw_include_all.h>

#include <fc2d_cudaclaw.h>

#include <fc2d_cudaclaw_cuda.h>
#include <fc2d_cudaclaw_options.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include <cudaclaw_user_fort.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* fix syntax highlighted */
#endif

typedef struct user_options
{
    double period;
    int claw_version;
    int cuda;
    int is_registered;

} user_options_t;

void swirl_link_solvers(fclaw_global_t *glob);

/* ------------------------------------- Options ---------------------------------------*/
user_options_t* swirl_options_register (fclaw_app_t * app,
                                        const char *configfile);

void swirl_options_store (fclaw_global_t* glob, user_options_t* user);

const user_options_t* swirl_get_options(fclaw_global_t* glob);


/* --------------------------------------- Cuda ----------------------------------------*/

void swirl_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);
void swirl_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2);
void swirl_assign_b4step2(cudaclaw_cuda_b4step2_t *b4step2);

void setprob();

/* ------------------------------------ Fortran ----------------------------------------*/
#if 0
#define SETPROB FCLAW_F77_FUNC(setprob, SETPROB)
void SETPROB();
#endif


#ifdef __cplusplus
}
#endif

#endif
