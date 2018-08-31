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

#ifndef CUDACLAW5_USER_FORT_H
#define CUDACLAW5_USER_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* --------------------------------------------------------------------------------
   Clawpack 5.0 routines

   These are provided for user convenience.  These files are not compiled
   into the library, but should be provided by the user.

   These signatures can be used if the user file matches these signatures 
   and subroutine name. Otherwise, the user should provide their own headers.
   ------------------------------------------------------------------------------- */

/* Macros for C/Fortran portability */

#define SETPROB            FCLAW_F77_FUNC(setprob,           SETPROB)
void SETPROB();


#define CUDACLAW5_SETPROB FCLAW_F77_FUNC(cudaclaw5_setprob, CUDACLAW5_SETPROB)
void CUDACLAW5_SETPROB();

#define CUDACLAW5_QINIT   FCLAW_F77_FUNC(cudaclaw5_qinit,   CUDACLAW5_QINIT)
void CUDACLAW5_QINIT(const int* meqn,const int* mbc,
                     const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy,
                     double q[], const int* maux, double aux[]);

#define CUDACLAW5_SETAUX  FCLAW_F77_FUNC(cudaclaw5_setaux,  CUDACLAW5_SETAUX)
void CUDACLAW5_SETAUX(const int* mbc,
                      const int* mx, const int* my,
                      const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      const int* maux, double aux[]);

#define CUDACLAW5_BC2     FCLAW_F77_FUNC(cudaclaw5_bc2,     CUDACLAW5_BC2)
void CUDACLAW5_BC2(const int* meqn, const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   const double q[], const int* maux,
                   const double aux[], const double* t,
                   const double* dt, const int mthbc[]);

#define CUDACLAW5_B4STEP2 FCLAW_F77_FUNC(cudaclaw5_b4step2, CUDACLAW5_B4STEP2)
void CUDACLAW5_B4STEP2(const int* mbc,
                       const int* mx, const int* my, const int* meqn,
                       double q[], const double* xlower,
                       const double* ylower,
                       const double* dx, const double* dy,
                       const double* t, const double* dt,
                       const int* maux, double aux[]);

#define CUDACLAW5_SRC2    FCLAW_F77_FUNC(cudaclaw5_src2,    CUDACLAW5_SRC2)
void CUDACLAW5_SRC2(const int* meqn,
                    const int* mbc, const int* mx,const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy, double q[],
                    const int* maux, double aux[], const double* t,
                    const double* dt);

/* Riemann solvers */
#define CUDACLAW5_RPN2    FCLAW_F77_FUNC(cudaclaw5_rpn2,    CUDACLAW5_RPN2)
void CUDACLAW5_RPN2(const int* ixy,const int* maxm, const int* meqn,
                    const int* mwaves, const int* maux,
                    const int* mbc,const int* mx,
                    double ql[], double qr[], double auxl[], double auxr[],
                    double wave[], double s[],double amdq[], double apdq[]);

#define CUDACLAW5_RPT2    FCLAW_F77_FUNC(cudaclaw5_rpt2,    CUDACLAW5_RPT2)
void CUDACLAW5_RPT2(const int* ixy, const int* imp,
                    const int* maxm, const int* meqn,
                    const int* mwaves, const int* maux,
                    const int* mbc,const int* mx,
                    double ql[], double qr[],
                    double aux1[], double aux2[],
                    double aux3[],  double asdq[],
                    double bmasdq[], double bpasdq[]);


#define CUDACLAW5_TAG4REFINEMENT FCLAW_F77_FUNC(cudaclaw5_tag4refinement, \
                                                          CUDACLAW5_TAG4REFINEMENT)

void CUDACLAW5_TAG4REFINEMENT(const int* mx,const int* my,
                              const int* mbc,const int* meqn,
                              const double* xlower, const double* ylower,
                              const double* dx, const double* dy,
                              const int* blockno,
                              double q[],
                              const double* tag_threshold,
                              const int* init_flag,
                              int* tag_patch);



#define CUDACLAW5_TAG4COARSENING FCLAW_F77_FUNC(cudaclaw5_tag4coarsening, \
                                                CUDACLAW5_TAG4COARSENING)

void CUDACLAW5_TAG4COARSENING(const int* mx, const int* my,
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
