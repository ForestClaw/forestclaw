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

#include <fclaw2d_include_all.h>

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

    int is_registered;

} user_options_t;

struct fclaw2d_global;
struct fclaw2d_patch;

#define SQUARE_SETPROB FCLAW_F77_FUNC(square_setprob, SQUARE_SETPROB)
void SQUARE_SETPROB();

void square_link_solvers(struct fclaw2d_global *glob);

/* ---------------------------------- Options ----------------------------------------- */

const user_options_t* square_get_options(struct fclaw2d_global* glob);

void square_options_store (fclaw2d_global_t* glob, user_options_t* user);


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

/* --------------------------------- Riemann Problems --------------------------------- */
void square_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* square_get_options(fclaw2d_global_t* glob);


/* ---------------------------- Fortran headers --------------------------------------- */

#define SQUARE_COMPUTE_ERROR FCLAW_F77_FUNC(square_compute_error,SQUARE_COMPUTE_ERROR)

void SQUARE_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define RPN2CONS_FW FCLAW_F77_FUNC(rpn2cons_fw, RPN2CONS_FW)
void RPN2CONS_FW(const int* ixy, const int* maxm, const int* meqn, 
                 const int* mwaves,const int* mbc, const int* mx, 
                 double ql[], double qr[],
                 double auxl[], double auxr[], double fwave[],
                 double s[], double amdq[], double apdq[]);

#define RPN2CONS_FW_MANIFOLD FCLAW_F77_FUNC(rpn2cons_fw_manifold, RPN2CONS_FW_MANIFOLD)
void RPN2CONS_FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                          const int* mwaves,
                          const int* mbc, const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double fwave[],
                          double s[], double amdq[], double apdq[]);


#define RPT2CONS FCLAW_F77_FUNC(rpt2cons, RPT2CONS)
void RPT2CONS(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
              const int* mbc, const int* mx, double ql[], double qr[],
              double aux1[], double aux2[], double aux3[], const int* imp,
              double dsdq[], double bmasdq[], double bpasdq[]);

#define RPT2CONS_MANIFOLD FCLAW_F77_FUNC(rpt2cons_manifold, RPT2CONS_MANIFOLD)
void RPT2CONS_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                       const int* mbc, const int* mx, double ql[], double qr[],
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


#define SQUARE_SETAUX FCLAW_F77_FUNC(square_setaux, SQUARE_SETAUX)

void SQUARE_SETAUX(const int* blockno, const int* mx, const int* my,
                   const int* mbc, const double* xlower, const double* ylower,
                   const double* dx, const double* dy, 
                   double area[],double edgelengths[],
                   double xp[], double yp[], double zp[],
                   double aux[],const int* maux);


#define SQUARE_SET_VELOCITIES FCLAW_F77_FUNC(square_set_velocities, \
                                             SQUARE_SET_VELOCITIES)

void SQUARE_SET_VELOCITIES(const int* blockno, const int* mx, const int* my,
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

#define  SQUARE_FORT_WRITE_FILE FCLAW_F77_FUNC(square_fort_write_file,  \
                                                SQUARE_FORT_WRITE_FILE)
void     SQUARE_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define SQUARE_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(square_fort_header_ascii, \
                        SQUARE_FORT_HEADER_ASCII)
void SQUARE_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);


#define SQUARE_TAG4REFINEMENT FCLAW_F77_FUNC(square_tag4refinement, \
                                              SQUARE_TAG4REFINEMENT)
void  SQUARE_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  SQUARE_TAG4COARSENING FCLAW_F77_FUNC(square_tag4coarsening, \
                                              SQUARE_TAG4COARSENING)
void  SQUARE_TAG4COARSENING(const int* mx, const int* my,
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


#define PERIODIC_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(periodic_fort_interpolate_face, \
                             PERIODIC_FORT_INTERPOLATE_FACE)
void PERIODIC_FORT_INTERPOLATE_FACE(const int* mx, const int* my, 
                                    const int* mbc,const int* meqn,
                                    double qcoarse[],double qfine[],
                                    const int* idir, const int* iside,
                                    const int* num_neighbors,
                                    const int* refratio, const int* igrid,
                                    struct fclaw2d_patch_transform_data** 
                                    transform_cptr);


#define PERIODIC_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(periodic_fort_interpolate_corner, \
                     PERIODIC_FORT_INTERPOLATE_CORNER)
void PERIODIC_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, 
                                      const int* mbc,const int* meqn, 
                                      const int* a_refratio, double this_q[],
                                      double neighbor_q[], const int* a_corner,
                                      struct fclaw2d_patch_transform_data** 
                                      transform_cptr);

#define PERIODIC_FORT_INTERPOLATE2FINE \
           FCLAW_F77_FUNC(periodic_fort_interpolate2fine, \
                          PERIODIC_FORT_INTERPOLATE2FINE)
void PERIODIC_FORT_INTERPOLATE2FINE(const int* mx,const int* my,
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
