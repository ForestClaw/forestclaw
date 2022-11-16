/*
Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef HEAT_USER_H
#define HEAT_USER_H

#include <fclaw2d_include_all.h>

#include <fc2d_thunderegg.h>

#include "heat_options.h"


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* --------------------------- Problem dependent functions -----------------------------*/

void heat_link_solvers(fclaw2d_global_t *glob);

void heat_run(fclaw2d_global_t *glob);


/* --------------------------- Fortran functions ---------------------------------------*/



#define HEAT_SETPROB FCLAW_F77_FUNC(heat_setprob,HEAT_SETPROB)

void HEAT_SETPROB();


#define HEAT_FORT_RHS FCLAW_F77_FUNC(heat_fort_rhs,HEAT_FORT_RHS)

void HEAT_FORT_RHS(const int* blockno, const int* mbc, const int* mx, 
                     const int* my, const int* meqn, const int* mfields, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, 
                     double *dt, int* method, double q[], double rhs[]);



#define HEAT_INIT FCLAW_F77_FUNC(heat_init,HEAT_INIT)

void HEAT_INIT(const int* blockno, const int* mbc, const int* mx, 
               const int* my, const int* meqn, 
               const double *xlower, const double *ylower,
               const double* dx, const double* dy, 
               double q[]);


#define HEAT_UPDATE_Q FCLAW_F77_FUNC(heat_update_q,HEAT_UPDATE_Q)

void HEAT_UPDATE_Q(const int* mbc, const int* mx, 
                   const int* my, const int* meqn, 
                   const int* mfields, double rhs[],
                   double q[]);



#define HEAT_FORT_BETA FCLAW_F77_FUNC(heat_fort_beta,HEAT_FORT_BETA)

void HEAT_FORT_BETA(const double* x, const double* y, const double* b, double grad[]);

#define HEAT_COMPUTE_ERROR FCLAW_F77_FUNC(heat_compute_error,HEAT_COMPUTE_ERROR)

void HEAT_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* mfields,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define HEAT_FORT_APPLY_BC FCLAW_F77_FUNC(heat_fort_apply_bc, \
                                            HEAT_FORT_APPLY_BC)

void HEAT_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                          const  int* mbc, const  int* mfields, 
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, const double* t,
                          int intersects_bc[], int mthbc[], 
                          double rhs[], fc2d_thunderegg_fort_eval_bc_t g_bc, 
                          int* cons_check, double flux_sum[]);


#define HEAT_FORT_EVAL_BC FCLAW_F77_FUNC(heat_fort_eval_bc, HEAT_FORT_EVAL_BC)

double HEAT_FORT_EVAL_BC(const int* iface, const double* t,const double* x, const double* y);


#define HEAT_NEUMANN FCLAW_F77_FUNC(heat_neumann, HEAT_NEUMANN)

double HEAT_NEUMANN(const int* iface, const double* t,const double* x, const double* y);

#define HEAT_TAG4REFINEMENT FCLAW_F77_FUNC(heat_tag4refinement,HEAT_TAG4REFINEMENT)

void HEAT_TAG4REFINEMENT(const int* mx,const int* my,
					const int* mbc,const int* meqn,
					const double* xlower, const double* ylower,
					const double* dx, const double* dy,
					const int* blockno,
					double q[],
					const double* tag_threshold,
					const int* init_flag,
					int* tag_patch);

#define HEAT_TAG4COARSENING FCLAW_F77_FUNC(heat_tag4coarsening,HEAT_TAG4COARSENING)

void HEAT_TAG4COARSENING(const int* mx, const int* my,
					const int* mbc, const int* meqn,
					double xlower[], double ylower[],
					const double* dx, const double* dy,
					const int* blockno,
					double q0[],double q1[],
					double q2[],double q3[],
					const double* tag_threshold,
                    const int* initflag,
					int* tag_patch);



/* ----------------------------- Fortran - output functions --------------------------- */

#define  HEAT_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(heat_fort_output_ascii, \
                          HEAT_FORT_OUTPUT_ASCII)
void HEAT_FORT_OUTPUT_ASCII(char* matname1,
                              int* mx,        int* my,
                              int* meqn,      int* mbc,
                              double* xlower, double* ylower,
                              double* dx,     double* dy,
                              double q[],double soln[], double error[],
                              int* patch_num, int* level,
                              int* blockno,   int* mpirank);

#define HEAT_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(heat_fort_header_ascii, \
                        HEAT_FORT_HEADER_ASCII)
void HEAT_FORT_HEADER_ASCII(char* matname1, char* matname2,
                              double* time, int* meqn, int* maux, 
                              int* ngrids);


#define HEAT_FORT_BC2 FCLAW_F77_FUNC(heat_fort_bc2, HEAT_FORT_BC2)

void HEAT_FORT_BC2(const int* meqn, const int* mbc, 
                   const int* mx, const int* my, 
                   const double *xlower, const double *ylower, 
                   const double *dx, const double *dy, 
                   double q[], double* t, double *dt, 
                   int intersects_bc[]);

fclaw2d_domain_t* heat_create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt);

void heat_run_program(fclaw2d_global_t* glob);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
