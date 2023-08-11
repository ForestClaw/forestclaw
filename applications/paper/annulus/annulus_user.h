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
    double revs_per_s;
    double cart_speed;
    double amplitude;
    double freq;
    double init_radius;
    double beta;    /* Ratio of inner radius to outer radius */

    const char *theta_string;
    double *theta;  /* Theta range */

    int refine_pattern;

    int claw_version;

    int is_registered;
}
user_options_t;


void annulus_link_solvers(fclaw_global_t *glob);

#define SETPROB_ANNULUS FCLAW_F77_FUNC(setprob_annulus, \
                                       SETPROB_ANNULUS)
void SETPROB_ANNULUS(void); 



user_options_t* annulus_options_register (fclaw_app_t * app,
                                          const char *configfile);

void annulus_options_store (fclaw_global_t* glob, user_options_t* user);

const user_options_t* annulus_get_options(fclaw_global_t *glob);

fclaw2d_map_context_t *
    fclaw2d_map_new_annulus (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double beta,
                             const double theta[]);


/* ----------------------
   Clawpack 4.6 headers
   ---------------------- */

#define ANNULUS_SETAUX FCLAW_F77_FUNC(annulus_setaux, ANNULUS_SETAUX)

void ANNULUS_SETAUX(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* xlower, const double* ylower,
                   const double* dx, const double* dy, 
                   double area[],double edgelengths[],
                   double xp[], double yp[], double zp[],
                   double aux[],const int* maux);


#define ANNULUS_SET_VELOCITIES FCLAW_F77_FUNC(annulus_set_velocities, \
                                             ANNULUS_SET_VELOCITIES)

void ANNULUS_SET_VELOCITIES(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* dx, const double* dy,
                   const double* xlower, const double* ylower,
                   const double *t, double xnormals[],double ynormals[],
                   double surfnormals[], double aux[],const int* maux);

#if 0
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
#endif                    


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


#define RPN2QAD_FLUX FCLAW_F77_FUNC(rpn2qad_flux, RPN2QAD_FLUX)
void RPN2QAD_FLUX(const int* meqn, const int* maux, const int* idir,
                               const int* iface, double q[], 
                               double auxvec[], double aux_edge[],
                               double flux[]);


#define RPN2QAD_ZERO FCLAW_F77_FUNC(rpn2qad_zero,RPN2QAD_ZERO)

void RPN2QAD_ZERO (const int* meqn, const int* maux, const int* idir,
                   const int* iface, double q[], double aux_center[], 
                   double aux_edge[], double flux[]);


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
