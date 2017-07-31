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
    int rp_solver;
    int example;
    int is_registered;

} user_options_t;

#define SWIRL_SETPROB FCLAW_F77_FUNC(swirl_setprob, SWIRL_SETPROB)
void SWIRL_SETPROB(int* example);

void swirlcons_link_solvers(fclaw2d_global_t *glob);

void swirlcons_problem_setup(fclaw2d_global_t* glob);

const user_options_t* swirlcons_get_options(fclaw2d_global_t* glob);

#define RPN2CONS_QS FCLAW_F77_FUNC(rpn2cons_qs,RPN2CONS_QS)
void RPN2CONS_QS(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
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

#define RPN2CONS_FW FCLAW_F77_FUNC(rpn2cons_fw, RPN2CONS_FW)
void RPN2CONS_FW(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                 const int* mbc, const int* mx, double ql[], double qr[],
                 double auxl[], double auxr[], double fwave[],
                 double s[], double amdq[], double apdq[]);

#define RPT2CONS FCLAW_F77_FUNC(rpt2cons, RPT2CONS)
void RPT2CONS(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
              const int* mbc, const int* mx, double ql[], double qr[],
              double aux1[], double aux2[], double aux3[], const int* imp,
              double dsdq[], double bmasdq[], double bpasdq[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
