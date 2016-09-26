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

#ifndef FC2D_SWIRL46_USER_FORT_H
#define FC2D_SWIRL46_USER_FORT_H

#include <fclaw2d_forestclaw.h>

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
#define SWIRL46_SETPROB FCLAW_F77_FUNC(swirl46_setprob,SWIRL46_SETPROB)
#define SWIRL46_QINIT   FCLAW_F77_FUNC(swirl46_qinit,  SWIRL46_QINIT)
#define SWIRL46_SETAUX  FCLAW_F77_FUNC(swirl46_setaux, SWIRL46_SETAUX)
#define SWIRL46_B4STEP2 FCLAW_F77_FUNC(swirl46_b4step2,SWIRL46_B4STEP2)
#define SWIRL46_SRC2    FCLAW_F77_FUNC(swirl46_src2,   SWIRL46_SRC2)
#define SWIRL46_BC2     FCLAW_F77_FUNC(swirl46_bc2,    SWIRL46_bc2)
#define SWIRL46_RPN2    FCLAW_F77_FUNC(swirl46_rpn2,   SWIRL46_RPN2)
#define SWIRL46_RPT2    FCLAW_F77_FUNC(swirl46_rpt2,   SWIRL46_RPT2)

void SWIRL46_SETPROB();

void SWIRL46_QINIT(const int* maxmx, const int* maxmy, const int* meqn,
                           const int* mbc, const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           double q[], const int* maux, double aux[]);

void SWIRL46_SETAUX(const int* maxmx, const int* maxmy, const int* mbc,
                            const int* mx, const int* my,
                            const double* xlower, const double* ylower,
                            const double* dx, const double* dy,
                            const int* maux, double aux[]);

void SWIRL46_BC2(const int* maxmx, const int* maxmy, const int* meqn,
                         const int* mbc, const int* mx, const int* my,
                         const double* xlower, const double* ylower,
                         const double* dx, const double* dy, const double q[],
                         const int* maux, const double aux[], const double* t,
                         const double* dt, const int mthbc[]);


void SWIRL46_B4STEP2(const int* maxmx, const int* maxmy, const int* mbc,
                             const int* mx, const int* my, const int* meqn,
                             double q[], const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const double* t, const double* dt,
                             const int* maux, double aux[]);

void SWIRL46_SRC2(const int* maxmx, const int* maxmy, const int* meqn,
                          const int* mbc, const int* mx,const int* my,
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, double q[],
                          const int* maux, double aux[], const double* t,
                          const double* dt);

void SWIRL46_RPN2(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                          const int* mbc,const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double wave[],
                          double s[], double amdq[], double apdq[]);

void SWIRL46_RPT2(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                          const int* mbc, const int* mx, double ql[], double qr[],
                          double aux1[], double aux2[], double aux3[], const int* imp,
                          double dsdq[], double bmasdq[], double bpasdq[]);


#define SWIRL46_TAG4REFINEMENT FCLAW_F77_FUNC(swirl46_tag4refinement,SWIRL46_TAG4REFINEMENT)
void SWIRL46_TAG4REFINEMENT(const int* mx,const int* my,
                    const int* mbc,const int* meqn,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* blockno,
                    double q[],
                    const double* tag_threshold,
                    const int* init_flag,
                    int* tag_patch);

#define SWIRL46_TAG4COARSENING FCLAW_F77_FUNC(swirl46_tag4coarsening,SWIRL46_TAG4COARSENING)
void SWIRL46_TAG4COARSENING(const int* mx, const int* my,
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
