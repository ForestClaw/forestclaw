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

#ifndef ANNULUS_USER_H
#define ANNULUS_USER_H

#include <fclaw2d_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
#endif    

typedef struct user_options
{
    int example;
    int mapping; 
    int initchoice;  /* Smooth or non-smooth */
    double revs_per_s;
    double twist;
    double vertical_speed;
    double init_radius;
    int color_equation;
    int use_stream;
    double beta;    /* Ratio of inner radius to outer radius */
    int refine_pattern;

    int claw_version;

    int is_registered;
}
user_options_t;


void annulus_link_solvers(fclaw2d_global_t *glob);

void annulus_problem_setup(fclaw2d_global_t *glob);

#define SETPROB_ANNULUS FCLAW_F77_FUNC(setprob_annulus, \
                                       SETPROB_ANNULUS)

void SETPROB_ANNULUS(void); 

#if 0
void SETPROB_ANNULUS(const int* example, 
                     const int* mapping, 
                     const int* initial_condition,
                     const double* revs_per_s,
                     const int *ceqn_in,
                     const int *use_stream_in,
                     const double* beta,
                     const int* refine_pattern);
#endif                     


void annulus_patch_setup(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx);

user_options_t* annulus_options_register (fclaw_app_t * app,
                                          const char *configfile);

void annulus_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* annulus_get_options(fclaw2d_global_t *glob);

fclaw2d_map_context_t *
    fclaw2d_map_new_annulus (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double beta,
                             const double twist,
                             const int mapping);


/* ----------------------
   Clawpack 4.6 headers
   ---------------------- */
#if 1
#define TORUS46_COMPUTE_ERROR FCLAW_F77_FUNC(torus46_compute_error,TORUS46_COMPUTE_ERROR)

void TORUS46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);
#endif

#define ANNULUS46_COMPUTE_ERROR FCLAW_F77_FUNC(annulus46_compute_error, \
                                             ANNULUS46_COMPUTE_ERROR)

void ANNULUS46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);



#define ANNULUS46_SETAUX  FCLAW_F77_FUNC(annulus46_setaux, ANNULUS46_SETAUX)
void ANNULUS46_SETAUX(const int* mbc, const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* maux, double aux[], const int* blockno,
                    double area[], 
                    double xd[], double yd[], double zd[],
                    double edgelengths[], 
                    double xnormals[], double ynormals[], 
                    double xtangents[], double ytangents[], 
                    double surfnormals[]);


#define  ANNULUS46_FORT_WRITE_FILE FCLAW_F77_FUNC(annulus46_fort_write_file,  \
                                                ANNULUS46_FORT_WRITE_FILE)
void     ANNULUS46_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define ANNULUS46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(annulus46_fort_header_ascii, \
                        ANNULUS46_FORT_HEADER_ASCII)
void ANNULUS46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);


#if 0
#define ANNULUS46_TAG4REFINEMENT FCLAW_F77_FUNC(annulus46_tag4refinement, \
                                              ANNULUS46_TAG4REFINEMENT)
void  ANNULUS46_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  ANNULUS46_TAG4COARSENING FCLAW_F77_FUNC(annulus46_tag4coarsening, \
                                              ANNULUS46_TAG4COARSENING)
void  ANNULUS46_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             int* tag_patch);
#endif



#define RPN2CONS_FW_MANIFOLD FCLAW_F77_FUNC(rpn2cons_fw_manifold, RPN2CONS_FW_MANIFOLD)
void RPN2CONS_FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                          const int* mwaves,
                          const int* mbc, const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double fwave[],
                          double s[], double amdq[], double apdq[]);


#define RPT2CONS_MANIFOLD FCLAW_F77_FUNC(rpt2cons_manifold, RPT2CONS_MANIFOLD)
void RPT2CONS_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                       const int* mbc, const int* mx, double ql[], double qr[],
                       double aux1[], double aux2[], double aux3[], const int* imp,
                       double dsdq[], double bmasdq[], double bpasdq[]);


#define RPN2_CONS_UPDATE_MANIFOLD FCLAW_F77_FUNC(rpn2_cons_update_manifold, \
                                                 RPN2_CONS_UPDATE_MANIFOLD)

void RPN2_CONS_UPDATE_MANIFOLD(const int* meqn, const int* maux, const int* idir,
                               const int* iface, double q[], 
                               double aux_center[], double aux_edge[],
                               double flux[]);

#define RPN2_CONS_UPDATE_ZERO FCLAW_F77_FUNC(rpn2_cons_update_zero, \
                                             RPN2_CONS_UPDATE_ZERO)

void RPN2_CONS_UPDATE_ZERO(const int* meqn, const int* maux, const int* idir,
                           const int* iface,
                           double q[], double aux_center[], double aux_edge[],
                           double flux[]);


#define ANNULUS46_RPT2ADV_MANIFOLD FCLAW_F77_FUNC(annulus46_rpt2adv_manifold, \
                                                   ANNULUS46_RPT2ADV_MANIFOLD)
void ANNULUS46_RPT2ADV_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                                 const int* mwaves,
                                 const int* mbc, const int* mx, double ql[], double qr[],
                                 double aux1[], double aux2[], double aux3[], const int* imp,
                                 double dsdq[], double bmasdq[], double bpasdq[]);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* ANNULUS_USER_H */
