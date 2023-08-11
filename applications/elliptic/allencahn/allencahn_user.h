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

#ifndef ALLENCAHN_USER_H
#define ALLENCAHN_USER_H

#include <fclaw2d_include_all.h>

#include <fc2d_thunderegg.h>

#include "allencahn_options.h"


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* --------------------------- Problem dependent functions -----------------------------*/

void allencahn_link_solvers(fclaw_global_t *glob);

void allencahn_run(fclaw_global_t *glob);


/* --------------------------- Fortran functions ---------------------------------------*/



#define ALLENCAHN_SETPROB FCLAW_F77_FUNC(allencahn_setprob,ALLENCAHN_SETPROB)

void ALLENCAHN_SETPROB();


#define ALLENCAHN_FORT_RHS FCLAW_F77_FUNC(allencahn_fort_rhs,ALLENCAHN_FORT_RHS)

void ALLENCAHN_FORT_RHS(const int* blockno, const int* mbc, const int* mx, 
                     const int* my, const int* meqn, const int* mfields, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, 
                     double *dt, int* method, double q[], double rhs[]);



#define ALLENCAHN_INIT FCLAW_F77_FUNC(allencahn_init,ALLENCAHN_INIT)

void ALLENCAHN_INIT(const int* meqn, const int* mbc, const int* mx, 
               const int* my,  
               const double *xlower, const double *ylower,
               const double* dx, const double* dy, 
               double q[]);


#define ALLENCAHN_UPDATE_Q FCLAW_F77_FUNC(allencahn_update_q,ALLENCAHN_UPDATE_Q)

void ALLENCAHN_UPDATE_Q(const int* mbc, const int* mx, 
                   const int* my, const int* meqn, 
                   const int* mfields, double rhs[],
                   double q[]);



#define ALLENCAHN_FORT_BETA FCLAW_F77_FUNC(allencahn_fort_beta,ALLENCAHN_FORT_BETA)

void ALLENCAHN_FORT_BETA(const double* x, const double* y, const double* b, double grad[]);

#define ALLENCAHN_COMPUTE_ERROR FCLAW_F77_FUNC(allencahn_compute_error,ALLENCAHN_COMPUTE_ERROR)

void ALLENCAHN_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* mfields,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define ALLENCAHN_FORT_APPLY_BC FCLAW_F77_FUNC(allencahn_fort_apply_bc, \
                                            ALLENCAHN_FORT_APPLY_BC)

void ALLENCAHN_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                          const  int* mbc, const  int* mfields, 
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, const double* t,
                          int intersects_bc[], int mthbc[], 
                          double rhs[], fc2d_thunderegg_fort_eval_bc_t g_bc, 
                          int* cons_check, double flux_sum[]);


#define ALLENCAHN_FORT_EVAL_BC FCLAW_F77_FUNC(allencahn_fort_eval_bc, ALLENCAHN_FORT_EVAL_BC)

double ALLENCAHN_FORT_EVAL_BC(const int* iface, const double* t,const double* x, const double* y);


#define ALLENCAHN_DIRICHLET FCLAW_F77_FUNC(allencahn_dirichlet, ALLENCAHN_DIRICHLET)

double ALLENCAHN_DIRICHLET(const int* iface, const double* t,const double* x, const double* y);



#define ALLENCAHN_FORT_BC2 FCLAW_F77_FUNC(allencahn_fort_bc2, ALLENCAHN_FORT_BC2)

void ALLENCAHN_FORT_BC2(const int* meqn, const int* mbc, 
                        const int* mx, const int* my, 
                        const double *xlower, const double *ylower, 
                        const double *dx, const double *dy, 
                        double q[], double* t, double *dt, 
                        int intersects_bc[]);


/* ----------------------------- Fortran - output functions --------------------------- */

#define  ALLENCAHN_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(allencahn_fort_output_ascii, \
                          ALLENCAHN_FORT_OUTPUT_ASCII)
void ALLENCAHN_FORT_OUTPUT_ASCII(char* matname1,
                              int* mx,        int* my,
                              int* meqn,      int* mbc,
                              double* xlower, double* ylower,
                              double* dx,     double* dy,
                              double q[],double soln[], double error[],
                              int* patch_num, int* level,
                              int* blockno,   int* mpirank);

#define ALLENCAHN_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(allencahn_fort_header_ascii, \
                        ALLENCAHN_FORT_HEADER_ASCII)
void ALLENCAHN_FORT_HEADER_ASCII(char* matname1, char* matname2,
                              double* time, int* meqn, int* maux, 
                              int* ngrids);




#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
