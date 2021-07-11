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

#ifndef SPHERE_USER_H
#define SPHERE_USER_H

#include "../all/transport_user.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int example;
    int mapping;
    int initial_condition;
    int refine_pattern;

    double *omega;
    const char* omega_string;

    int claw_version;

    int is_registered;

} user_options_t;

struct fclaw2d_global;
struct fclaw2d_patch;

#if 0
/* So syntax highlighting works */
#endif

#if 0
#define SPHERE_SETPROB FCLAW_F77_FUNC(sphere_setprob, SPHERE_SETPROB)
void SPHERE_SETPROB();
#endif

void sphere_link_solvers(struct fclaw2d_global *glob);

/* ---------------------------------- Options ----------------------------------------- */

const user_options_t* sphere_get_options(struct fclaw2d_global* glob);

void sphere_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* sphere_options_register (fclaw_app_t * app, const char *configfile);

/* --------------------------------- Mappings ----------------------------------------- */
#if 0
fclaw2d_map_context_t * fclaw2d_map_new_cubedsphere(const double scale[], 
                                                    const double rotate[]);

fclaw2d_map_context_t * fclaw2d_map_new_pillowsphere(const double scale[]);
#endif

/* ----------------------------------- Options ---------------------------------------- */

void sphere_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* sphere_get_options(fclaw2d_global_t* glob);

/* ---------------------------------- Compute Error ----------------------------------- */

#define SPHERE46_COMPUTE_ERROR FCLAW_F77_FUNC(sphere46_compute_error, \
                                              SPHERE46_COMPUTE_ERROR)

void SPHERE46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);

#define SPHERE5_COMPUTE_ERROR FCLAW_F77_FUNC(sphere5_compute_error, \
                                              SPHERE5_COMPUTE_ERROR)

void SPHERE5_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);

/* ---------------------------------- Riemann Solvers --------------------------------- */

#define RPN2CONS_FW_MANIFOLD FCLAW_F77_FUNC(rpn2cons_fw_manifold, RPN2CONS_FW_MANIFOLD)
void RPN2CONS_FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                          const int* mwaves, 
                          const int* mbc, const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double fwave[],
                          double s[], double amdq[], double apdq[], const int* maux);


#define RPT2CONS_MANIFOLD FCLAW_F77_FUNC(rpt2cons_manifold, RPT2CONS_MANIFOLD)
void RPT2CONS_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                       const int* maux, const int* mbc, const int* mx, 
                       double ql[], double qr[],
                       double aux1[], double aux2[], double aux3[], const int* imp,
                       double dsdq[], double bmasdq[], double bpasdq[]);


#define RPN2_CONS_UPDATE FCLAW_F77_FUNC(rpn2_cons_update,RPN2_CONS_UPDATE)

void RPN2_CONS_UPDATE(const int* meqn, const int* maux, const int* idir, const int* iface,
                      double q[], double aux_center[], double aux_edge[], double flux[]);


#define RPN2_CONS_UPDATE_MANIFOLD FCLAW_F77_FUNC(rpn2_cons_update_manifold, \
                                                 RPN2_CONS_UPDATE_MANIFOLD)

void RPN2_CONS_UPDATE_MANIFOLD(const int* meqn, const int* maux, const int* idir,
                               const int* iface,
                               double q[], double aux_center[], double aux_edge[],
                               double flux[]);

/* -------------------------------------- Output -------------------------------------- */

#define  SPHERE46_FORT_WRITE_FILE FCLAW_F77_FUNC(sphere46_fort_write_file,  \
                                                SPHERE46_FORT_WRITE_FILE)
void  SPHERE46_FORT_WRITE_FILE(char* matname1,
                             int* mx,        int* my,
                             int* meqn,      int* mbc,
                             double* xlower, double* ylower,
                             double* dx,     double* dy,
                             double q[],     double error[], double soln[],
                             double *time,
                             int* patch_num, int* level,
                             int* blockno,   int* mpirank);

#define SPHERE46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(sphere46_fort_header_ascii, \
                        SPHERE46_FORT_HEADER_ASCII)
void SPHERE46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#define  SPHERE5_FORT_WRITE_FILE FCLAW_F77_FUNC(sphere5_fort_write_file,  \
                                                SPHERE5_FORT_WRITE_FILE)
void  SPHERE5_FORT_WRITE_FILE(char* matname1,
                             int* mx,        int* my,
                             int* meqn,      int* mbc,
                             double* xlower, double* ylower,
                             double* dx,     double* dy,
                             double q[],     double error[], double soln[],
                             double *time,
                             int* patch_num, int* level,
                             int* blockno,   int* mpirank);

#define SPHERE5_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(sphere5_fort_header_ascii, \
                        SPHERE5_FORT_HEADER_ASCII)
void SPHERE5_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

/* -------------------------------------- Old -------------------------------------- */

#if 0
#define SPHERE_BASIS_COMPLETE FCLAW_F77_FUNC(sphere_basis_complete, SPHERE_BASIS_COMPLETE)

void SPHERE_BASIS_COMPLETE(const double* x, const double *y,
                           double t[], double tinv[], double uderivs[], 
                           const int* flag);
#endif                           


#define SPHERE_SETAUX FCLAW_F77_FUNC(sphere_setaux, SPHERE_SETAUX)


#if 0
void SPHERE_SETAUX(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* xlower, const double* ylower,
                   const double* dx, const double* dy, 
                   double area[],double edgelengths[],
                   double xp[], double yp[], double zp[],
                   double aux[],const int* maux);


#define SPHERE_SET_VELOCITIES FCLAW_F77_FUNC(sphere_set_velocities, \
                                             SPHERE_SET_VELOCITIES)

void SPHERE_SET_VELOCITIES(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* dx, const double* dy,
                   const double* xlower, const double* ylower,
                   const double *t, double xnormals[],double ynormals[],
                   double surfnormals[], double aux[],const int* maux);
#endif

#define  SPHERE_FORT_WRITE_FILE FCLAW_F77_FUNC(sphere_fort_write_file,  \
                                                SPHERE_FORT_WRITE_FILE)
void  SPHERE_FORT_WRITE_FILE(char* matname1,
                             int* mx,        int* my,
                             int* meqn,      int* mbc,
                             double* xlower, double* ylower,
                             double* dx,     double* dy,
                             double q[],     double error[], double soln[],
                             double *time,
                             int* patch_num, int* level,
                             int* blockno,   int* mpirank);

#define SPHERE_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(sphere_fort_header_ascii, \
                        SPHERE_FORT_HEADER_ASCII)
void SPHERE_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#define  SPHERE_FORT_WRITE_FILE FCLAW_F77_FUNC(sphere_fort_write_file,  \
                                                SPHERE_FORT_WRITE_FILE)
void  SPHERE_FORT_WRITE_FILE(char* matname1,
                             int* mx,        int* my,
                             int* meqn,      int* mbc,
                             double* xlower, double* ylower,
                             double* dx,     double* dy,
                             double q[],     double error[], double soln[],
                             double *time,
                             int* patch_num, int* level,
                             int* blockno,   int* mpirank);

#define SPHERE_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(sphere_fort_header_ascii, \
                        SPHERE_FORT_HEADER_ASCII)
void SPHERE_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#if 0
#define SPHERE_TAG4REFINEMENT FCLAW_F77_FUNC(sphere_tag4refinement, \
                                              SPHERE_TAG4REFINEMENT)
void  SPHERE_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const double *t, const int* blockno,
                             double q[], const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  SPHERE_TAG4COARSENING FCLAW_F77_FUNC(sphere_tag4coarsening, \
                                              SPHERE_TAG4COARSENING)
void  SPHERE_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const double *t, 
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             int* tag_patch);

#endif

#define USER_EXCEEDS_THRESHOLD FCLAW_F77_FUNC(user_exceeds_threshold, \
                                              USER_EXCEEDS_THRESHOLD)

int USER_EXCEEDS_THRESHOLD(int* blockno,
                           double qval[], 
                           double* qmin, double *qmax,
                           double quad[], 
                           double *dx, double *dy, 
                           double *xc, double *yc, 
                           int* tag_threshold, 
                           int* init_flag,
                           int* is_ghost);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
