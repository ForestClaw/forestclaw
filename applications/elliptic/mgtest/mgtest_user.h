/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef MGTEST_USER_H
#define MGTEST_USER_H

#include <fclaw2d_include_all.h>
#include "mgtest_options.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* --------------------------- Problem dependent functions -----------------------------*/

void mgtest_link_solvers(fclaw2d_global_t *glob);

void mgtest_problem_setup(fclaw2d_global_t* glob);

/* --------------------------- Fortran functions ---------------------------------------*/

#define MGTEST_FORT_RHS FCLAW_F77_FUNC(mgtest_fort_rhs,MGTEST_FORT_RHS)

void MGTEST_FORT_RHS(const int* mbc, const int* mx, const int* my, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, double q[]);


#define MGTEST_SETPROB FCLAW_F77_FUNC(mgtest_setprob,MGTEST_SETPROB)

void MGTEST_SETPROB(const int* rhs_choice, const double *alpha,
                    const double* x0, const double* y0,
                    const double* a,  const double* b);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
