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

#ifndef PERIODIC_USER_H
#define PERIODIC_USER_H

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
    double uvel;
    double vvel;
    int claw_version;
    int is_registered;

} user_options_t;

void periodic_link_solvers(fclaw2d_global_t *glob);

void periodic_problem_setup(fclaw2d_global_t* glob);

/* ------------------------------------- Options ---------------------------------------*/
user_options_t* periodic_options_register (fclaw_app_t * app,
                                        const char *configfile);

void periodic_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* periodic_get_options(fclaw2d_global_t* glob);

/* ------------------------------------ Fortran ----------------------------------------*/
#define PERIODIC_SETPROB FCLAW_F77_FUNC(periodic_setprob, PERIODIC_SETPROB)
void PERIODIC_SETPROB();

#define CLAWPACK46_TAG4REFINEMENT FCLAW_F77_FUNC(clawpack46_tag4refinement, \
                                                 CLAWPACK46_TAG4REFINEMENT)
void CLAWPACK46_TAG4REFINEMENT(const int* mx,const int* my,
                               const int* mbc,const int* meqn,
                               const double* xlower, const double* ylower,
                               const double* dx, const double* dy,
                               const int* blockno,
                               double q[],
                               const double* tag_threshold,
                               const int* init_flag,
                               int* tag_patch);

#define CLAWPACK5_TAG4COARSENING FCLAW_F77_FUNC(clawpack5_tag4coarsening, \
                                                CLAWPACK5_TAG4COARSENING)
void CLAWPACK5_TAG4COARSENING(const int* mx, const int* my,
                              const int* mbc, const int* meqn,
                              const double* xlower, const double* ylower,
                              const double* dx, const double* dy,
                              const int* blockno,
                              double q0[],double q1[],
                              double q2[],double q3[],
                              const double* tag_threshold,
                              const int* initflag,
                              int* tag_patch);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
