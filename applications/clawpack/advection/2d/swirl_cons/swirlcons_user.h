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

#ifndef SWIRLCONS_USER_H
#define SWIRLCONS_USER_H

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
    double period;
    int rp_solver;
    int example;
    int is_registered;

} user_options_t;

#define SWIRL_SETPROB FCLAW_F77_FUNC(swirl_setprob, SWIRL_SETPROB)
void SWIRL_SETPROB(double* tperiod, int* example);

void swirlcons_link_solvers(fclaw2d_global_t *glob);

void swirlcons_problem_setup(fclaw2d_global_t* glob);

const user_options_t* swirlcons_get_options(fclaw2d_global_t* glob);

/* Mappings */
fclaw2d_map_context_t* fclaw2d_map_new_nomap();

#define RPN2CONS_QSTAR FCLAW_F77_FUNC(rpn2cons_qstar,RPN2CONS_QSTAR)
void RPN2CONS_QSTAR(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                    const int* mbc,const int* mx, double ql[], double qr[],
                    double auxl[], double auxr[], double wave[],
                    double s[], double amdq[], double apdq[]);

#define RPN2CONS_WD FCLAW_F77_FUNC(rpn2cons_wd,RPN2CONS_WD)
void RPN2CONS_WD(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
              const int* mbc,const int* mx, double ql[], double qr[],
              double auxl[], double auxr[], double wave[],
              double s[], double amdq[], double apdq[]);

#define RPN2CONS_EC FCLAW_F77_FUNC(rpn2cons_ec,RPN2CONS_EC)
void RPN2CONS_EC(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                 const int* mbc,const int* mx, double ql[], double qr[],
                 double auxl[], double auxr[], double wave[],
                 double s[], double amdq[], double apdq[]);

#define RPN2FWAVE FCLAW_F77_FUNC(rpn2fwave, RPN2FWAVE)
void RPN2FWAVE(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
               const int* mbc, const int* mx, double ql[], double qr[],
               double auxl[], double auxr[], double fwave[],
               double s[], double amdq[], double apdq[]);

#define RPT2CONS_CC FCLAW_F77_FUNC(rpt2cons_cc, RPT2CONS_CC)
void RPT2CONS_CC(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                 const int* mbc, const int* mx, double ql[], double qr[],
                 double aux1[], double aux2[], double aux3[], const int* imp,
                 double dsdq[], double bmasdq[], double bpasdq[]);

#define SWIRLCONS_BC2 FCLAW_F77_FUNC(swirlcons_bc2,SWIRLCONS_BC2)
void SWIRLCONS_BC2(const int* maxmx, const int* maxmy, const int* meqn,
                     const int* mbc, const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy, const double q[],
                     const int* maux, const double aux[], const double* t,
                     const double* dt, const int mthbc[]);

#define SWIRLCONS_B4STEP2 FCLAW_F77_FUNC(swirlcons_b4step2,SWIRLCONS_B4STEP2)
void SWIRLCONS_B4STEP2(const int* maxmx, const int* maxmy, const int* mbc,
                       const int* mx, const int* my, const int* meqn,
                       double q[], const double* xlower, const double* ylower,
                       const double* dx, const double* dy,
                       const double* t, const double* dt,
                       const int* maux, double aux[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
