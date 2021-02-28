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

#ifndef SGN_FORT_H
#define SGN_FORT_H

#include <fclaw2d_include_all.h>
#include <fc2d_thunderegg.h>


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
    /* Fix syntax highlighting */
#endif    


#if 0
#define TSUNAMI_SETPROB FCLAW_F77_FUNC(tsunami_setprob, TSUNAMI_SETPROB)
void TSUNAMI_SETPROB();
#endif


/* ***************************** FORTRAN - Riemann solvers **************************** */


#define SGN_FORT_RHS FCLAW_F77_FUNC(sgn_fort_rhs,SGN_FORT_RHS)

void SGN_FORT_RHS(const int* blockno, const int* mbc, const int* mx, 
                     const int* my, const int* meqn, const int* mfields, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, 
                     double q[], double rhs[]);


#define SGN_UPDATE_Q FCLAW_F77_FUNC(sgn_update_q,SGN_UPDATE_Q)

void SGN_UPDATE_Q(const int* mx, const int *my, const int* mbc, const int* meqn,
                  const int* mfields, const double* xlower, const double* ylower,
                  const double* dx, const double* dy, 
                  const double *t, const double* dt, 
                  const int* maux, 
                  double aux[], double q[], double D[]);


#define SGN_FORT_APPLY_BC FCLAW_F77_FUNC(sgn_fort_apply_bc, \
                                            SGN_FORT_APPLY_BC)

void SGN_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                          const  int* mbc, const  int* mfields, 
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, const double* t,
                          int intersects_bc[], int mthbc[], 
                          double rhs[], fc2d_thunderegg_fort_eval_bc_t g_bc, 
                          int* cons_check, double flux_sum[]);


#define SGN_FORT_EVAL_BC FCLAW_F77_FUNC(sgn_fort_eval_bc, SGN_FORT_EVAL_BC)

double SGN_FORT_EVAL_BC(const int* iface, const double* t,const double* x, const double* y);


#define SGN_NEUMANN FCLAW_F77_FUNC(sgn_neumann, SGN_NEUMANN)

double SGN_NEUMANN(const int* iface, const double* t,const double* x, const double* y);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
