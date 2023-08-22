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

#ifndef CLAWPACK5_USER_FORT_H
#define CLAWPACK5_USER_FORT_H

#include <fclaw_forestclaw.h>
#include <fclaw_package.h>

#include "fc3d_clawpack5_options.h"

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
   compiled into the library, but will be provided by the user.
   -------------------------------------------------------------------- */

/* Macros for C/Fortran portability */
#define CLAWPACK5_SETPROB FCLAW_F77_FUNC(clawpack5_setprob, CLAWPACK5_SETPROB)
#define CLAWPACK5_QINIT   FCLAW_F77_FUNC(clawpack5_qinit,   CLAWPACK5_QINIT)
#define CLAWPACK5_SETAUX  FCLAW_F77_FUNC(clawpack5_setaux,  CLAWPACK5_SETAUX)
#define CLAWPACK5_B4STEP2 FCLAW_F77_FUNC(clawpack5_b4step2, CLAWPACK5_B4STEP2)
#define CLAWPACK5_SRC2    FCLAW_F77_FUNC(clawpack5_src2,    CLAWPACK5_SRC2)
#define CLAWPACK5_BC2     FCLAW_F77_FUNC(clawpack5_bc2,     CLAWPACK5_BC2)
#define CLAWPACK5_RPN2    FCLAW_F77_FUNC(clawpack5_rpn2,    CLAWPACK5_RPN2)
#define CLAWPACK5_RPT2    FCLAW_F77_FUNC(clawpack5_rpt2,    CLAWPACK5_RPT2)

/* These will be converted to MACROS slowly ... */

void CLAWPACK5_SETPROB();

void CLAWPACK5_QINIT(const int* meqn,const int* mbc,
                     const int* mx, const int* my, const int* mz,
                     const double* xlower, const double* ylower, const double* zlower,
                     const double* dx, const double* dy, const double* dz,
                     double q[], const int* maux, double aux[]);

void CLAWPACK5_SETAUX(const int* mbc,
                      const int* mx, const int* my,
                      const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      const int* maux, double aux[]);

void CLAWPACK5_BC2(const int* meqn, const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   const double q[], const int* maux,
                   const double aux[], const double* t,
                   const double* dt, const int mthbc[]);

void CLAWPACK5_B4STEP2(const int* mbc,
                       const int* mx, const int* my, const int* meqn,
                       double q[], const double* xlower,
                       const double* ylower,
                       const double* dx, const double* dy,
                       const double* t, const double* dt,
                       const int* maux, double aux[]);

void CLAWPACK5_SRC2(const int* meqn,
                    const int* mbc, const int* mx,const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy, double q[],
                    const int* maux, double aux[], const double* t,
                    const double* dt);

/* Riemann solvers */
void CLAWPACK5_RPN2(const int* ixy,const int* maxm, const int* meqn,
                    const int* mwaves, const int* maux,
                    const int* mbc,const int* mx,
                    double ql[], double qr[], double auxl[], double auxr[],
                    double wave[], double s[],double amdq[], double apdq[]);

void CLAWPACK5_RPT2(const int* ixy, const int* imp,
                    const int* maxm, const int* meqn,
                    const int* mwaves, const int* maux,
                    const int* mbc,const int* mx,
                    double ql[], double qr[],
                    double aux1[], double aux2[],
                    double aux3[],  double asdq[],
                    double bmasdq[], double bpasdq[]);


#define CLAWPACK5_TAG4REFINEMENT FCLAW_F77_FUNC(clawpack5_tag4refinement, \
                                                          CLAWPACK5_TAG4REFINEMENT)

void CLAWPACK5_TAG4REFINEMENT(const int* mx,const int* my,
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
                              int* tag_patch);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
