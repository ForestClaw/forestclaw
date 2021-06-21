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

#ifndef ANNULUS_USER_H
#define ANNULUS_USER_H

#include <fclaw2d_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef struct user_options
{
    int example;
    double beta;  /* Ratio of inner radius to outer radius */

    const char *theta_string;
    double *theta;  /* Theta range */

    int claw_version;

    int is_registered;
}
user_options_t;


void annulus_link_solvers(fclaw2d_global_t *glob);

void annulus_problem_setup(fclaw2d_global_t *glob);

#define SETPROB_ANNULUS FCLAW_F77_FUNC(setprob_annulus,SETPROB_ANNULUS)
void SETPROB_ANNULUS();

void annulus_patch_setup(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx);

user_options_t* annulus_options_register (fclaw_app_t * app,
                                          const char *configfile);

void annulus_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* annulus_get_options(fclaw2d_global_t *glob);

fclaw2d_map_context_t *
    fclaw2d_map_new_annulus (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double alpha,
                             const double theta[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* ANNULUS_USER_H */
