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

#ifndef SWIRLCONS_USER_H
#define SWIRLCONS_USER_H

#include <fclaw2d_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int rp_solver;
    int example;
    int mapping;

    double alpha;    /* Used by five-patch mapping */

    int initial_condition;  /* Smooth or non-smooth */
    int color_equation;
    int use_stream;

    double *center;
    const char *center_string;

    int is_registered;

} user_options_t;

struct fclaw2d_global;
struct fclaw2d_patch;

#define SWIRLCONS_SETPROB FCLAW_F77_FUNC(swirlcons_setprob, SWIRLCONS_SETPROB)
void SWIRLCONS_SETPROB();

void swirlcons_link_solvers(struct fclaw2d_global *glob);

void swirlcons_problem_setup(struct fclaw2d_global* glob);

/* ---------------------------------- Options ----------------------------------------- */

const user_options_t* swirlcons_get_options(struct fclaw2d_global* glob);

void swirlcons_options_store (fclaw2d_global_t* glob, user_options_t* user);


user_options_t* swirlcons_options_register (fclaw_app_t * app,
                                       const char *configfile);

/* --------------------------------- Mappings ----------------------------------------- */
fclaw_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double rotate[],
                                                 const double alpha);

fclaw_map_context_t* fclaw2d_map_new_cart (fclaw_map_context_t* brick,
                                             const double scale[],
                                             const double shift[],
                                             const double rotate[]);

fclaw_map_context_t* fclaw2d_map_new_bilinear(fclaw_map_context_t *brick,
                                                const double scale[],
                                                const double shift[],
                                                const double rotate[],
                                                const double center[]);

void swirlcons_patch_setup_manifold(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch *this_patch,
                                    int this_block_idx,
                                    int this_patch_idx);

/* --------------------------------- Riemann Problems --------------------------------- */
void swirlcons_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* swirlcons_get_options(fclaw2d_global_t* glob);


/* ---------------------------- Fortran headers --------------------------------------- */

#define SWIRL46_COMPUTE_ERROR FCLAW_F77_FUNC(swirl46_compute_error,SWIRL46_COMPUTE_ERROR)

void SWIRL46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define RPN2CONS_QS FCLAW_F77_FUNC(rpn2cons_qs,RPN2CONS_QS)
void RPN2CONS_QS(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                 const int* mbc,const int* mx, double ql[], double qr[],
                 double auxl[], double auxr[], double wave[],
                 double s[], double amdq[], double apdq[]);

#define RPN2CONS_QS_MANIFOLD FCLAW_F77_FUNC(rpn2cons_qs_manifold,RPN2CONS_QS_MANIFOLD)
void RPN2CONS_QS_MANIFOLD(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                          const int* mbc,const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double wave[],
                          double s[], double amdq[], double apdq[]);

#define RPN2CONS_WD FCLAW_F77_FUNC(rpn2cons_wd,RPN2CONS_WD)
void RPN2CONS_WD(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
              const int* mbc,const int* mx, double ql[], double qr[],
              double auxl[], double auxr[], double wave[],
              double s[], double amdq[], double apdq[]);

#define RPN2CONS_EC FCLAW_F77_FUNC(rpn2cons_ec,RPN2CONS_EC)
void RPN2CONS_EC(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                 const int* mbc,const int* mx, double ql[], double qr[],
                 double auxl[], double auxr[], double wave[],
                 double s[], double amdq[], double apdq[]);

#define RPN2CONS_EC_MANIFOLD FCLAW_F77_FUNC(rpn2cons_ec_manifold,RPN2CONS_EC_MANIFOLD)
void RPN2CONS_EC_MANIFOLD(const int* ixy,const int* maxm, const int* meqn, 
                          const int* mwaves,const int* mbc,const int* mx, 
                          double ql[], double qr[],
                          double auxl[], double auxr[], double wave[],
                          double s[], double amdq[], double apdq[]);


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


#define CLAWPACK46_SETAUX_MANIFOLD FCLAW_F77_FUNC(clawpack46_setaux_manifold, \
                                               CLAWPACK46_SETAUX_MANIFOLD)

void CLAWPACK46_SETAUX_MANIFOLD(const int* mbc,
                            const int* mx, const int* my,
                            const double* xlower, const double* ylower,
                            const double* dx, const double* dy,
                            const int* maux, double aux[],
                            const int* blockno,
                            double xp[], double yp[], double zp[],
                            double area[],double edgelengths[],
                            double xnormals[],double ynormals[],
                            double surfnormals[]);

#define SWIRL46_SETAUX FCLAW_F77_FUNC(swirl46_setaux, SWIRL46_SETAUX)

void SWIRL46_SETAUX(const int* mbc,
                            const int* mx, const int* my,
                            const double* xlower, const double* ylower,
                            const double* dx, const double* dy,
                            const int* maux, double aux[],
                            const int* blockno,
                            double area[],double edgelengths[],
                            double xnormals[],double ynormals[],
                            double surfnormals[]);


#define  SWIRL46_FORT_WRITE_FILE FCLAW_F77_FUNC(swirl46_fort_write_file,  \
                                                SWIRL46_FORT_WRITE_FILE)
void     SWIRL46_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define SWIRL46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(swirl46_fort_header_ascii, \
                        SWIRL46_FORT_HEADER_ASCII)
void SWIRL46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);


#define SWIRL46_TAG4REFINEMENT FCLAW_F77_FUNC(swirl46_tag4refinement, \
                                              SWIRL46_TAG4REFINEMENT)
void  SWIRL46_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  SWIRL46_TAG4COARSENING FCLAW_F77_FUNC(swirl46_tag4coarsening, \
                                              SWIRL46_TAG4COARSENING)
void  SWIRL46_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             int* tag_patch);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
