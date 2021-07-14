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

#ifndef CLAWPACK46_USER_FORT_H
#define CLAWPACK46_USER_FORT_H

#ifdef __cplusplus
extern "C"
{
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

/* ------------------------------- User defined routines ------------------------------ */

#define CLAWPACK46_SETPROB FCLAW_F77_FUNC(clawpack_setprob,CLAWPACK46_SETPROB)
void CLAWPACK46_SETPROB();

#define CLAWPACK46_QINIT   FCLAW_F77_FUNC(clawpack46_qinit,  CLAWPACK46_QINIT)
void CLAWPACK46_QINIT(const int* meqn, const int* mbc,
                      const int* mx, const int* my, const int* mz,
                      const double* xlower, const double* ylower,
                      const double* zlower, 
                      const double* dx, const double* dy, const double *dz,
                      double q[], const int* maux, double aux[]);
  
#define CLAWPACK46_SETAUX  FCLAW_F77_FUNC(clawpack46_setaux, CLAWPACK46_SETAUX)
void CLAWPACK46_SETAUX(const int* mbc, const int* mx, const int* my, const int* mz, 
                       const double* xlower, const double* ylower, const double* zlower,
                       const double* dx, const double* dy, const double* dz,
                       const int* maux, double aux[]);


#define CLAWPACK46_BC3     FCLAW_F77_FUNC(clawpack46_bc3,    CLAWPACK46_bc3)
void CLAWPACK46_BC3(const int* meqn, const int* mbc, 
                    const int* mx, const int* my, const int* mz,  
                    const double* xlower, const double* ylower, const double* zlower,
                    const double* dx, const double* dy, const double *dz,
                    const double q[], const int* maux, const double aux[], 
                    const double* t, const double* dt, const int mthbc[]);


#define CLAWPACK46_B4STEP3 FCLAW_F77_FUNC(clawpack46_b4step3,CLAWPACK46_B4STEP3)
void CLAWPACK46_B4STEP3(const int* mbc, const int* mx, const int* my, const int* mz,
                        const int* meqn,
                        double q[], const double* xlower, const double* ylower,
                        const double* zlower,
                        const double* dx, const double* dy, const double* dz,
                        const double* t, const double* dt,
                        const int* maux, double aux[]);

#define CLAWPACK46_SRC3    FCLAW_F77_FUNC(clawpack46_src3,   CLAWPACK46_SRC3)
void CLAWPACK46_SRC3(const int* meqn,
                     const int* mbc, const int* mx,const int* my, const int* mz,
                     const double* xlower, const double* ylower, const double* zlower,
                     const double* dx, const double* dy, const double* dz,
                     double q[],
                     const int* maux, double aux[], const double* t,
                     const double* dt);

/* ---------------------------------- Riemann solvers  -------------------------------- */

#define CLAWPACK46_RPN3 FCLAW_F77_FUNC(clawpack46_rpn3,  CLAWPACK46_RPN3)
void CLAWPACK46_RPN3(const int* ixyz,const int* maxm, 
                     const int* meqn, const int* mwaves,
                     const int* maux, const int* mbc, const int* mx, 
                     double ql[], double qr[],
                     double auxl[], double auxr[], double wave[],
                     double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT3 FCLAW_F77_FUNC(clawpack46_rpt3, CLAWPACK46_RPT3)
void CLAWPACK46_RPT3(const int* ixyz, const int* icoor, const int* imp,
                     const int *maxm, const int* meqn, const int* mwaves, 
                     const int *maux, 
                     const int* mbc, const int* mx, double ql[], double qr[],
                     double aux1[], double aux2[], double aux3[], 
                     double asdq[], double bmasdq[], double bpasdq[]);

#define CLAWPACK46_RPTT3    FCLAW_F77_FUNC(clawpack46_rptt3, CLAWPACK46_RPTT3)
void CLAWPACK46_RPTT3(const int* ixyz, const int* icoor, const int* imp,
                      const int* impt, const int* maxm, const int* meqn,
                      const int* mwaves, const int* maux,
                      const int* mbc,const int* mx,
                      double ql[], double qr[],
                      double aux1[], double aux2[],
                      double aux3[],  double bsasdq[],
                      double cmbsasdq[], double cpbsasdq[]);

/* --------------------------------- Tagging routines  -------------------------------- */


#define CLAWPACK46_TAG4REFINEMENT FCLAW_F77_FUNC(clawpack46_tag4refinement, \
                                                 CLAWPACK46_TAG4REFINEMENT)

void CLAWPACK46_TAG4REFINEMENT(const int* mx,const int* my, const int* mz,
                               const int* mbc,const int* meqn,
                               const double* xlower, const double* ylower,
                               const double* zlower, 
                               const double* dx, const double* dy, const int* dz,
                               const int* blockno, double q[],
                               const double* tag_threshold,
                               const int* init_flag,
                               int* tag_patch);


#define CLAWPACK46_TAG4COARSENING FCLAW_F77_FUNC(clawpack46_tag4coarsening, \
                                                CLAWPACK46_TAG4COARSENING)

void CLAWPACK46_TAG4COARSENING(const int* mx, const int* my, const *mz,
                               const int* mbc, const int* meqn,
                               const double* xlower, const double* ylower,
                               const int* zlower, 
                               const double* dx, const double* dy, const int* dz,
                               const int* blockno,
                               double *qfine[], const double* tag_threshold,
                               int* tag_patch);

#ifdef __cplusplus
}
#endif

#endif
