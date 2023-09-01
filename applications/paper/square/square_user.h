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

#ifndef SQUARE_USER_H
#define SQUARE_USER_H

#include <fclaw_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
/* Needed so syntax highlighing works */
#endif


typedef struct user_options
{
    int example;
    int mapping;

    int initial_condition;  /* Smooth or non-smooth */

    double alpha;    /* Used by five-patch mapping */

    double *velocity;
    const char *velocity_string;

    double *center;
    const char *center_string;

    int claw_version;

    int is_registered;

} user_options_t;

struct fclaw_global;
struct fclaw_patch_t;

#define SQUARE_SETPROB FCLAW_F77_FUNC(square_setprob, SQUARE_SETPROB)
void SQUARE_SETPROB();

void square_link_solvers(struct fclaw_global *glob);

/* ---------------------------------- Options ----------------------------------------- */

const user_options_t* square_get_options(struct fclaw_global* glob);

void square_options_store (fclaw_global_t* glob, user_options_t* user);


user_options_t* square_options_register (fclaw_app_t * app,
                                       const char *configfile);

/* --------------------------------- Mappings ----------------------------------------- */
fclaw2d_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double alpha);

fclaw2d_map_context_t* fclaw2d_map_new_cart (fclaw2d_map_context_t* brick,
                                             const double scale[],
                                             const double shift[]);

fclaw2d_map_context_t* fclaw2d_map_new_bilinear(fclaw2d_map_context_t *brick,
                                                const double scale[],
                                                const double shift[],
                                                const double center[]);

fclaw2d_map_context_t* fclaw2d_map_new_identity(fclaw2d_map_context_t *brick);

/* ---------------------------------- Compute Error ------------------------------------- */

#define SQUARE46_COMPUTE_ERROR FCLAW_F77_FUNC(square46_compute_error,SQUARE46_COMPUTE_ERROR)

void SQUARE46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define SQUARE5_COMPUTE_ERROR FCLAW_F77_FUNC(square5_compute_error,SQUARE5_COMPUTE_ERROR)

void SQUARE5_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);

/* ------------------------- clawpack 4.6 - Riemann solvers ---------------------------- */

#define CLAWPACK46_RPN2FW_MANIFOLD FCLAW_F77_FUNC(clawpack46_rpn2fw_manifold, \
                                                  CLAWPACK46_RPN2FW_MANIFOLD)
void CLAWPACK46_RPN2FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                                const int* mwaves, const int* mbc, const int* mx, 
                                double ql[], double qr[],
                                double auxl[], double auxr[], double fwave[],
                                double s[], double amdq[], double apdq[],
                                const int* maux);

/* Note : maux is passed in as argument 5 */
#define CLAWPACK46_RPT2FW_MANIFOLD FCLAW_F77_FUNC(clawpack46_rpt2fw_manifold, \
                                                CLAWPACK46_RPT2fw_MANIFOLD)
void CLAWPACK46_RPT2FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                              const int* mwaves, const int* maux, 
                              const int* mbc, const int* mx, double ql[], double qr[],
                              double aux1[], double aux2[], double aux3[], 
                              const int* imp, double dsdq[], 
                              double bmasdq[], double bpasdq[]);

/* ------------------------- clawpack 5.0 - Riemann solvers ---------------------------- */

/* Both flux2 and flux2 using fwaves are in a single function */

#define CLAWPACK5_RPN2FW_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpn2fw_manifold, \
                                                 CLAWPACK5_RPN2FW_MANIFOLD)
void CLAWPACK5_RPN2FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                                const int* mwaves, const int* maux, 
                                const int* mbc, const int* mx, 
                                double ql[], double qr[],
                                double auxl[], double auxr[], double fwave[],
                                double s[], double amdq[], double apdq[]);


#define CLAWPACK5_RPT2_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpt2_manifold, \
                                               CLAWPACK5_RPT2_MANIFOLD)

void CLAWPACK5_RPT2_MANIFOLD(const int* ixy, const int* imp, const int* maxm, 
                             const int* meqn, const int* mwaves, const int* maux,
                             const int* mbc, const int* mx, double ql[], double qr[],
                             double aux1[], double aux2[], double aux3[], 
                             double dsdq[], double bmasdq[], double bpasdq[]);


/* ------------------------------ General - Riemann solvers --------------------------------- */

#define RPN2_CONS_UPDATE FCLAW_F77_FUNC(rpn2_cons_update,RPN2_CONS_UPDATE)

void RPN2_CONS_UPDATE(const int* meqn, const int* maux, const int* idir, const int* iface,
                      double q[], double aux_center[], double aux_edge[], double flux[]);


#define RPN2_CONS_UPDATE_MANIFOLD FCLAW_F77_FUNC(rpn2_cons_update_manifold, \
                                                 RPN2_CONS_UPDATE_MANIFOLD)

void RPN2_CONS_UPDATE_MANIFOLD(const int* meqn, const int* maux, const int* idir,
                               const int* iface,
                               double q[], double aux_center[], double aux_edge[],
                               double flux[]);

#define RPN2_CONS_UPDATE_ZERO FCLAW_F77_FUNC(rpn2_cons_update_zero, \
                                                 RPN2_CONS_UPDATE_ZERO)

void RPN2_CONS_UPDATE_ZERO(const int* meqn, const int* maux, const int* idir,
                               const int* iface,
                               double q[], double aux_center[], double aux_edge[],
                               double flux[]);

/* ---------------------------------------- SETAUX ------------------------------------------ */
#define SQUARE46_SETAUX FCLAW_F77_FUNC(square46_setaux, SQUARE46_SETAUX)

void SQUARE46_SETAUX(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* xlower, const double* ylower,
                   const double* dx, const double* dy, 
                   double area[],double edgelengths[],
                   double xp[], double yp[], double zp[],
                   double aux[],const int* maux);


#define SQUARE5_SETAUX FCLAW_F77_FUNC(square5_setaux, SQUARE5_SETAUX)

void SQUARE5_SETAUX(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* xlower, const double* ylower,
                   const double* dx, const double* dy, 
                   double area[],double edgelengths[],
                   double xp[], double yp[], double zp[],
                   double aux[],const int* maux);

/* ---------------------------------------- SET_VELOCITIES --------------------------------------- */

#define SQUARE46_SET_VELOCITIES FCLAW_F77_FUNC(square46_set_velocities, \
                                             SQUARE46_SET_VELOCITIES)

void SQUARE46_SET_VELOCITIES(const int* blockno, const int* mx, const int* my,
                             const int* mbc, const double* dx, const double* dy,
                             const double* xlower, const double* ylower,
                             const double *t, double xnormals[],double ynormals[],
                             double surfnormals[], double aux[],const int* maux);


#define SQUARE5_SET_VELOCITIES FCLAW_F77_FUNC(square5_set_velocities, \
                                             SQUARE5_SET_VELOCITIES)

void SQUARE5_SET_VELOCITIES(const int* blockno, const int* mx, const int* my,
                            const int* mbc, const double* dx, const double* dy,
                            const double* xlower, const double* ylower,
                            const double *t, double xnormals[],double ynormals[],
                            double surfnormals[], double aux[],const int* maux);


#if 0
#define CART_BASIS_COMPLETE FCLAW_F77_FUNC(cart_basis_complete, CART_BASIS_COMPLETE)

void CART_BASIS_COMPLETE(const int* blockno, const double* x, const double *y,
                          double t[], double tinv[], double uderivs[], 
                          const int* flag);


#define SQUARE_SETAUX FCLAW_F77_FUNC(square_setaux, SQUARE_SETAUX)

void SQUARE_SETAUX(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   double area[],double edgelengths[],
                   double xnormals[],double ynormals[],
                   double surfnormals[],
                   double aux[],const int* maux);
#endif

/* ---------------------------------------- OUTPUT --------------------------------------- */

#define  SQUARE46_FORT_WRITE_FILE FCLAW_F77_FUNC(square46_fort_write_file,  \
                                                SQUARE46_FORT_WRITE_FILE)
void     SQUARE46_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define SQUARE46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(square46_fort_header_ascii, \
                        SQUARE46_FORT_HEADER_ASCII)
void SQUARE46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#define  SQUARE5_FORT_WRITE_FILE FCLAW_F77_FUNC(square5_fort_write_file,  \
                                                SQUARE5_FORT_WRITE_FILE)
void     SQUARE5_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define SQUARE5_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(square5_fort_header_ascii, \
                        SQUARE5_FORT_HEADER_ASCII)
void SQUARE5_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);



/* ---------------------------------------- TAGGING --------------------------------------- */

#define SQUARE46_TAG4REFINEMENT FCLAW_F77_FUNC(square46_tag4refinement, \
                                              SQUARE46_TAG4REFINEMENT)
void  SQUARE46_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  SQUARE46_TAG4COARSENING FCLAW_F77_FUNC(square46_tag4coarsening, \
                                              SQUARE46_TAG4COARSENING)
void  SQUARE46_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             const int* initflag,
                             int* tag_patch);


#define SQUARE5_TAG4REFINEMENT FCLAW_F77_FUNC(square5_tag4refinement, \
                                              SQUARE5_TAG4REFINEMENT)
void  SQUARE5_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  SQUARE5_TAG4COARSENING FCLAW_F77_FUNC(square5_tag4coarsening, \
                                              SQUARE5_TAG4COARSENING)
void  SQUARE5_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             const int* initflag,
                             int* tag_patch);




#if 0
#define SQUARE_TAG4REFINEMENT_WAVELET FCLAW_F77_FUNC(square_tag4refinement_wavelet, \
                                              SQUARE_TAG4REFINEMENT_WAVELET)
void  SQUARE_TAG4REFINEMENT_WAVELET(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  SQUARE_TAG4COARSENING_WAVELET FCLAW_F77_FUNC(square_tag4coarsening_wavelet, \
                                              SQUARE_TAG4COARSENING_WAVELET)
void  SQUARE_TAG4COARSENING_WAVELET(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             int* tag_patch);


#define SQUARE_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(square_fort_interpolate_face, \
                             PERIODIC_FORT_INTERPOLATE_FACE)
void SQUARE_FORT_INTERPOLATE_FACE(const int* mx, const int* my, 
                                    const int* mbc,const int* meqn,
                                    double qcoarse[],double qfine[],
                                    const int* idir, const int* iside,
                                    const int* num_neighbors,
                                    const int* refratio, const int* igrid,
                                    struct fclaw2d_patch_transform_data** 
                                    transform_cptr);


#define SQUARE_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(square_fort_interpolate_corner, \
                     SQUARE_FORT_INTERPOLATE_CORNER)
void SQUARE_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, 
                                      const int* mbc,const int* meqn, 
                                      const int* a_refratio, double this_q[],
                                      double neighbor_q[], const int* a_corner,
                                      struct fclaw2d_patch_transform_data** 
                                      transform_cptr);

#define SQUARE_FORT_INTERPOLATE2FINE \
           FCLAW_F77_FUNC(square_fort_interpolate2fine, \
                          SQUARE_FORT_INTERPOLATE2FINE)
void SQUARE_FORT_INTERPOLATE2FINE(const int* mx,const int* my,
                                    const int* mbc, const int* meqn,
                                    double qcoarse[], double qfine[],
                                    double areacoarse[], double areafine[],
                                    const int* igrid, const int* manifold);
#endif                                    


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
