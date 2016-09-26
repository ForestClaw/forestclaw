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

#ifndef CLAWPACK46_USER_FORT_H
#define CLAWPACK46_USER_FORT_H

#include <fclaw2d_forestclaw.h>
#include <fclaw_package.h>

#include "fc2d_clawpack46_options.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* --------------------------------------------------------------------
   Classic routines
   - These are provided only for convenience;  these files are not
   compiled into the library, but should be provided by the user.

   Users maybe define these files how they wish. These signatures can
   be used if the user file matches these signatures and subroutine name.
   Otherwise, the user should provide their own headers.
   -------------------------------------------------------------------- */

/* Macros for C/Fortran portability */
#define SETPROB                 FCLAW_F77_FUNC(setprob,SETPROB)
#define CLAWPACK46_SETPROB FCLAW_F77_FUNC(clawpack46_setprob,CLAWPACK46_SETPROB)
#define CLAWPACK46_QINIT   FCLAW_F77_FUNC(clawpack46_qinit,  CLAWPACK46_QINIT)
#define CLAWPACK46_SETAUX  FCLAW_F77_FUNC(clawpack46_setaux, CLAWPACK46_SETAUX)
#define CLAWPACK46_B4STEP2 FCLAW_F77_FUNC(clawpack46_b4step2,CLAWPACK46_B4STEP2)
#define CLAWPACK46_SRC2    FCLAW_F77_FUNC(clawpack46_src2,   CLAWPACK46_SRC2)
#define CLAWPACK46_BC2     FCLAW_F77_FUNC(clawpack46_bc2,    CLAWPACK46_bc2)
#define CLAWPACK46_RPN2    FCLAW_F77_FUNC(clawpack46_rpn2,   CLAWPACK46_RPN2)
#define CLAWPACK46_RPT2    FCLAW_F77_FUNC(clawpack46_rpt2,   CLAWPACK46_RPT2)


void SETPROB();

void CLAWPACK46_SETPROB();

void CLAWPACK46_QINIT(const int* maxmx, const int* maxmy, const int* meqn,
                           const int* mbc, const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           double q[], const int* maux, double aux[]);

void CLAWPACK46_SETAUX(const int* maxmx, const int* maxmy, const int* mbc,
                            const int* mx, const int* my,
                            const double* xlower, const double* ylower,
                            const double* dx, const double* dy,
                            const int* maux, double aux[]);

void CLAWPACK46_BC2(const int* maxmx, const int* maxmy, const int* meqn,
                         const int* mbc, const int* mx, const int* my,
                         const double* xlower, const double* ylower,
                         const double* dx, const double* dy, const double q[],
                         const int* maux, const double aux[], const double* t,
                         const double* dt, const int mthbc[]);


void CLAWPACK46_B4STEP2(const int* maxmx, const int* maxmy, const int* mbc,
                             const int* mx, const int* my, const int* meqn,
                             double q[], const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const double* t, const double* dt,
                             const int* maux, double aux[]);

void CLAWPACK46_SRC2(const int* maxmx, const int* maxmy, const int* meqn,
                          const int* mbc, const int* mx,const int* my,
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, double q[],
                          const int* maux, double aux[], const double* t,
                          const double* dt);

void CLAWPACK46_RPN2(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                          const int* mbc,const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double wave[],
                          double s[], double amdq[], double apdq[]);

void CLAWPACK46_RPT2(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
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
