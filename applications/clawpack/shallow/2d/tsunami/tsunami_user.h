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

#ifndef TSUNAMI_USER_H
#define TSUNAMI_USER_H

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
    double g;
    double a;
    double b;
    double h0;

    double breaking;
    double alpha;
    double dry_tolerance;
    double sea_level;

    int claw_version;

    int is_registered;
} user_options_t;


#define TSUNAMI_SETPROB FCLAW_F77_FUNC(tsunami_setprob, TSUNAMI_SETPROB)
void TSUNAMI_SETPROB();


void tsunami_problem_setup(fclaw2d_global_t *glob);
void tsunami_link_solvers(fclaw2d_global_t *glob);

user_options_t* tsunami_options_register (fclaw_app_t * app,
                                          const char *configfile);

void tsunami_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* tsunami_get_options(fclaw2d_global_t* glob);

fclaw2d_map_context_t* fclaw2d_map_new_nomap();


/* ***************************** FORTRAN - Riemann solvers **************************** */

#define RPN2_TSUNAMI    FCLAW_F77_FUNC(rpn2_tsunami, RPN2_TSUNAMI)
void RPN2_TSUNAMI(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc,const int* mx, double ql[], double qr[],
                  double auxl[], double auxr[], double wave[],
                  double s[], double amdq[], double apdq[]);

#define RPT2_TSUNAMI    FCLAW_F77_FUNC(rpt2_tsunami, RPT2_TSUNAMI)
void RPT2_TSUNAMI(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc, const int* mx, double ql[], double qr[],
                  double aux1[], double aux2[], double aux3[], const int* imp,
                  double dsdq[], double bmasdq[], double bpasdq[]);

#define RPN2_GEOCLAW    FCLAW_F77_FUNC(rpn2_geoclaw,   RPN2_GEOCLAW)
void RPN2_GEOCLAW(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc,const int* mx, double ql[], double qr[],
                  double auxl[], double auxr[], double wave[],
                  double s[], double amdq[], double apdq[]);

#define RPT2_GEOCLAW    FCLAW_F77_FUNC(rpt2_geoclaw,   RPT2_GEOCLAW)
void RPT2_GEOCLAW(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
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
