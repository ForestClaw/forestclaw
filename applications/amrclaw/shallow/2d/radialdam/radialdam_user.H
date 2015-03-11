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

#ifndef RADIALDAM_USER_H
#define RADIALDAM_USER_H

#include "fclaw2d_forestclaw.H"
#include "fclaw2d_clawpatch.H"
#include "fc2d_clawpack46.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    double g;
    double x0;
    double y0;
    double r0;
    double hin;
    double hout;

    int is_registered;
} user_options_t;


#define RADIALDAM_SETPROB FCLAW_F77_FUNC(radialdam_setprob, RADIALDAM_SETPROB)
void RADIALDAM_SETPROB(const double *gamma, const double* x0, const double* y0,
                       const double* r0, const double* hin,
                       const double* hinf);

void radialdam_problem_setup(fclaw2d_domain_t* domain);
void radialdam_link_solvers(fclaw2d_domain_t *domain);


#define RPN2SW FCLAW_F77_FUNC(rpn2sw,RPN2SW)
void RPN2SW(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
            const int* mbc,const int* mx, double ql[], double qr[],
            double auxl[], double auxr[], double wave[],
            double s[], double amdq[], double apdq[]);

#define RPT2SW FCLAW_F77_FUNC(rpt2sw,RPT2SW)
void RPT2SW(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
            const int* mbc, const int* mx, double ql[], double qr[],
            double aux1[], double aux2[], double aux3[], const int* imp,
            double dsdq[], double bmasdq[], double bpasdq[]);

fclaw2d_map_context_t* fclaw2d_map_new_nomap();

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
