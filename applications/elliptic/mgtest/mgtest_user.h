/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef MGTEST_USER_H
#define MGTEST_USER_H

#include <fclaw2d_include_all.h>

#include <fc2d_multigrid.h>

#include "mgtest_options.h"


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* --------------------------- Problem dependent functions -----------------------------*/

void mgtest_link_solvers(fclaw2d_global_t *glob);

/* --------------------------- Fortran functions ---------------------------------------*/

#define MGTEST_SETPROB FCLAW_F77_FUNC(mgtest_setprob,MGTEST_SETPROB)

void MGTEST_SETPROB();


#define MGTEST_FORT_RHS FCLAW_F77_FUNC(mgtest_fort_rhs,MGTEST_FORT_RHS)

void MGTEST_FORT_RHS(const int* blockno, const int* mbc, const int* mx, 
                     const int* my, const int* mfields, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, double rhs[]);


#define MGTEST_FORT_BETA FCLAW_F77_FUNC(mgtest_beta,MGTEST_FORT_BETA)

void MGTEST_FORT_BETA(const double* x, const double* y, const double* b, double grad[]);

#define MGTEST_COMPUTE_ERROR FCLAW_F77_FUNC(mgtest_compute_error,MGTEST_COMPUTE_ERROR)

void MGTEST_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* mfields,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define MGTEST_FORT_APPLY_BC FCLAW_F77_FUNC(mgtest_fort_apply_bc, \
                                            MGTEST_FORT_APPLY_BC)

void MGTEST_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                          const  int* mbc, const  int* mfields, 
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, const double* t,
                          int intersects_bc[], int mthbc[], 
                          double rhs[], fc2d_multigrid_fort_eval_bc_t g_bc, 
                          int* cons_check, double flux_sum[]);


#define MGTEST_FORT_EVAL_BC FCLAW_F77_FUNC(mgtest_fort_eval_bc, MGTEST_FORT_EVAL_BC)

double MGTEST_FORT_EVAL_BC(const int* iface, const double* t,const double* x, const double* y);



/* ----------------------------- Fortran - output functions --------------------------- */

#define  MGTEST_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(mgtest_fort_output_ascii, \
                          MGTEST_FORT_OUTPUT_ASCII)
void MGTEST_FORT_OUTPUT_ASCII(char* matname1,
                              int* mx,        int* my,
                              int* meqn,      int* mbc,
                              double* xlower, double* ylower,
                              double* dx,     double* dy,
                              double q[],double soln[], double error[],
                              int* patch_num, int* level,
                              int* blockno,   int* mpirank);

#define MGTEST_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(mgtest_fort_header_ascii, \
                        MGTEST_FORT_HEADER_ASCII)
void MGTEST_FORT_HEADER_ASCII(char* matname1, char* matname2,
                              double* time, int* meqn, int* maux, 
                              int* ngrids);




#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
