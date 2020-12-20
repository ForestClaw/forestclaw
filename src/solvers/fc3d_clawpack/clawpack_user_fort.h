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

#ifndef CLAWPACK_USER_FORT_H
#define CLAWPACK_USER_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* --------------------------------------------------------------------------------
   Clawpack 4.6 routines

   These are provided for user convenience.  These files are not compiled
   into the library, but should be provided by the user.

   These signatures can be used if the user file matches these signatures 
   and subroutine name. Otherwise, the user should provide their own headers.
   ------------------------------------------------------------------------------- */



#define SETPROB            FCLAW_F77_FUNC(setprob,           SETPROB)
void SETPROB();

#define CLAWPACK_SETPROB FCLAW_F77_FUNC(clawpack_setprob,CLAWPACK_SETPROB)
void CLAWPACK_SETPROB();

#define CLAWPACK_QINIT   FCLAW_F77_FUNC(clawpack_qinit,  CLAWPACK_QINIT)
void CLAWPACK_QINIT(const int* maxmx, const int* maxmy, const int* meqn,
                      const int* mbc, const int* mx, const int* my,
                      const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      double q[], const int* maux, double aux[]);
  
#define CLAWPACK_SETAUX  FCLAW_F77_FUNC(clawpack_setaux, CLAWPACK_SETAUX)
void CLAWPACK_SETAUX(const int* maxmx, const int* maxmy, const int* mbc,
                       const int* mx, const int* my,
                       const double* xlower, const double* ylower,
                       const double* dx, const double* dy,
                       const int* maux, double aux[]);

#define CLAWPACK_BC2     FCLAW_F77_FUNC(clawpack_bc2,    CLAWPACK_bc2)
void CLAWPACK_BC2(const int* maxmx, const int* maxmy, const int* meqn,
                    const int* mbc, const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy, const double q[],
                    const int* maux, const double aux[], const double* t,
                    const double* dt, const int mthbc[]);


#define CLAWPACK_B4STEP2 FCLAW_F77_FUNC(clawpack_b4step2,CLAWPACK_B4STEP2)
void CLAWPACK_B4STEP2(const int* maxmx, const int* maxmy, const int* mbc,
                        const int* mx, const int* my, const int* meqn,
                        double q[], const double* xlower, const double* ylower,
                        const double* dx, const double* dy,
                        const double* t, const double* dt,
                        const int* maux, double aux[]);

#define CLAWPACK_SRC2    FCLAW_F77_FUNC(clawpack_src2,   CLAWPACK_SRC2)
void CLAWPACK_SRC2(const int* maxmx, const int* maxmy, const int* meqn,
                     const int* mbc, const int* mx,const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy, double q[],
                     const int* maux, double aux[], const double* t,
                     const double* dt);

#define CLAWPACK_RPN2    FCLAW_F77_FUNC(clawpack_rpn2,   CLAWPACK_RPN2)
void CLAWPACK_RPN2(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc,const int* mx, double ql[], double qr[],
                     double auxl[], double auxr[], double wave[],
                     double s[], double amdq[], double apdq[]);

#define CLAWPACK_RPT2    FCLAW_F77_FUNC(clawpack_rpt2,   CLAWPACK_RPT2)
void CLAWPACK_RPT2(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc, const int* mx, double ql[], double qr[],
                     double aux1[], double aux2[], double aux3[], const int* imp,
                     double dsdq[], double bmasdq[], double bpasdq[]);



#define CLAWPACK_TAG4REFINEMENT FCLAW_F77_FUNC(clawpack_tag4refinement, \
                                                 CLAWPACK_TAG4REFINEMENT)

void CLAWPACK_TAG4REFINEMENT(const int* mx,const int* my,
                               const int* mbc,const int* meqn,
                               const double* xlower, const double* ylower,
                               const double* dx, const double* dy,
                               const int* blockno,
                               double q[],
                               const double* tag_threshold,
                               const int* init_flag,
                               int* tag_patch);



#define CLAWPACK_TAG4COARSENING FCLAW_F77_FUNC(clawpack_tag4coarsening, \
                                                CLAWPACK_TAG4COARSENING)

void CLAWPACK_TAG4COARSENING(const int* mx, const int* my,
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
