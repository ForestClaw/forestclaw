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

#ifndef FC2D_CLAWPACK46_H
#define FC2D_CLAWPACK46_H

#include <forestclaw2d.h>
#include <fclaw_options.h>
#include <clawpack_fort.H>
#include <amr_forestclaw.H>

#include "fc2d_clawpack46_options.h"
#include <fclaw_package.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif



typedef void (*fc2d_clawpack46_setprob_t)();

typedef void (*fc2d_clawpack46_bc2_t)(const int* maxmx, const int* maxmy,
                                      const int* meqn, const int* mbc,
                                      const int* mx, const int* my,
                                      const double* xlower, const double* ylower,
                                      const double* dx, const double* dy,
                                      const double q[], const int* maux,
                                      const double aux[], const double* t,
                                      const double* dt, const int mthbc[]);

typedef  void (*fc2d_clawpack46_qinit_t)(const int* maxmx, const int* maxmy,
                                         const int* meqn,const int* mbc,
                                         const int* mx, const int* my,
                                         const double* xlower, const double* ylower,
                                         const double* dx, const double* dy,
                                         double q[], const int* maux, double aux[]);

typedef void (*fc2d_clawpack46_setaux_t)(const int* maxmx, const int* maxmy, const int* mbc,
                                         const int* mx, const int* my,
                                         const double* xlower, const double* ylower,
                                         const double* dx, const double* dy,
                                         const int* maux, double aux[]);

typedef void (*fc2d_clawpack46_b4step2_t)(const int* maxmx, const int* maxmy,
                                          const int* mbc,
                                          const int* mx, const int* my, const int* meqn,
                                          double q[], const double* xlower,
                                          const double* ylower,
                                          const double* dx, const double* dy,
                                          const double* t, const double* dt,
                                          const int* maux, double aux[]);

typedef void (*fc2d_clawpack46_src2_t)(const int* maxmx, const int* maxmy, const int* meqn,
                                       const int* mbc, const int* mx,const int* my,
                                       const double* xlower, const double* ylower,
                                       const double* dx, const double* dy, double q[],
                                       const int* maux, double aux[], const double* t,
                                       const double* dt);

typedef void (*fc2d_clawpack46_rpn2_t)(const int* ixy,const int* maxm, const int* meqn,
                                       const int* mwaves, const int* mbc,const int* mx,
                                       double ql[], double qr[], double auxl[], double auxr[],
                                       double wave[], double s[],double amdq[], double apdq[]);


typedef void (*fc2d_clawpack46_rpt2_t)(const int* ixy, const int* maxm, const int* meqn,
                                       const int* mwaves, const int* mbc,const int* mx,
                                       double ql[], double qr[], double aux1[], double aux2[],
                                       double aux3[], const int* imp, double dsdq[],
                                       double bmasdq[], double bpasdq[]);


typedef void (*fc2d_clawpack46_flux2_t)(const int* ixy,const int* maxm, const int* meqn,
                                        const int* maux,const int* mbc,const int* mx,
                                        double q1d[], double dtdx1d[],
                                        double aux1[], double aux2[], double aux3[],
                                        double faddm[],double faddp[], double gaddm[],
                                        double gaddp[],double cfl1d[], double fwave[],
                                        double s[], double amdq[],double apdq[],double cqxx[],
                                        double bmasdq[], double bpasdq[],
                                        fc2d_clawpack46_rpn2_t rpn2,
                                        fc2d_clawpack46_rpt2_t rpt2,
                                        const int* mwaves, const int* mcapa,
                                        int method[], int mthlim[]);

typedef struct fc2d_clawpack46_vtable
{
    fc2d_clawpack46_setprob_t setprob;
    fc2d_clawpack46_bc2_t bc2;
    fc2d_clawpack46_qinit_t qinit;
    fc2d_clawpack46_setaux_t setaux;
    fc2d_clawpack46_b4step2_t b4step2;
    fc2d_clawpack46_src2_t src2;
    fc2d_clawpack46_rpn2_t rpn2;
    fc2d_clawpack46_rpt2_t rpt2;
} fc2d_clawpack46_vtable_t;

void fc2d_clawpack46_set_vtable(const fc2d_clawpack46_vtable_t* vt);

#define CLAWPACK46_BC2 FCLAW_F77_FUNC(clawpack46_bc2,CLAWPACK46_BC2)
void CLAWPACK46_BC2(const int* maxmx, const int* maxmy, const int* meqn,
                     const int* mbc, const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy, const double q[],
                     const int* maux, const double aux[], const double* t,
                     const double* dt, const int mthbc[]);

/* --------------------------------------------------------------------
   Classic routines
   - These are provided only for convenience;  these files are not
   compiled into the library, but will be provided by the user.
   -------------------------------------------------------------------- */

/* Macros for C/Fortran portability */
#define SETPROB FCLAW_F77_FUNC(setprob,SETPROB)
#define QINIT   FCLAW_F77_FUNC(qinit,QINIT)
#define SETAUX  FCLAW_F77_FUNC(setaux,SETAUX)
#define B4STEP2 FCLAW_F77_FUNC(b4step2,B4STEP2)
#define SRC2    FCLAW_F77_FUNC(src2,SRC2)
#define BC2     FCLAW_F77_FUNC(bc2,BC2)
#define RPN2    FCLAW_F77_FUNC(rpn2,RPN2)
#define RPT2    FCLAW_F77_FUNC(rpt2,RPT2)

/* These will be converted to MACROS slowly ... */

void setprob_();

void qinit_(const int* maxmx, const int* maxmy, const int* meqn,
            const int* mbc, const int* mx, const int* my,
            const double* xlower, const double* ylower,
            const double* dx, const double* dy,
            double q[], const int* maux, double aux[]);

void setaux_(const int* maxmx, const int* maxmy, const int* mbc,
             const int* mx, const int* my,
             const double* xlower, const double* ylower,
             const double* dx, const double* dy,
             const int* maux, double aux[]);

void bc2_(const int* maxmx, const int* maxmy, const int* meqn,
          const int* mbc, const int* mx, const int* my,
          const double* xlower, const double* ylower,
          const double* dx, const double* dy, const double q[],
          const int* maux, const double aux[], const double* t,
          const double* dt, const int mthbc[]);


void b4step2_(const int* maxmx, const int* maxmy, const int* mbc,
              const int* mx, const int* my, const int* meqn,
              double q[], const double* xlower, const double* ylower,
              const double* dx, const double* dy,
              const double* t, const double* dt,
              const int* maux, double aux[]);

void src2_(const int* maxmx, const int* maxmy, const int* meqn,
           const int* mbc, const int* mx,const int* my,
           const double* xlower, const double* ylower,
           const double* dx, const double* dy, double q[],
           const int* maux, double aux[], const double* t,
           const double* dt);

/* Riemann solvers */
void rpn2_(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
           const int* mbc,const int* mx, double ql[], double qr[],
           double auxl[], double auxr[], double wave[],
           double s[], double amdq[], double apdq[]);

void rpt2_(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
           const int* mbc, const int* mx, double ql[], double qr[],
           double aux1[], double aux2[], double aux3[], const int* imp,
           double dsdq[], double bmasdq[], double bpasdq[]);

/* --------------------------------------------------------------------
   Time stepping
   -------------------------------------------------------------------- */
#define CLAWPACK46_STEP2_WRAP FCLAW_F77_FUNC(clawpack46_step2_wrap,CLAWPACK46_STEP2_WRAP)
void CLAWPACK46_STEP2_WRAP(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int method[], const int mthlim[],
                            const int* mcapa, const int* mwaves, const int* mx,
                            const int* my, double qold[], double auxold[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cfl, double work[], const int* mwork,
                            const double* xlower, const double* ylower, const int* level,
                            const double* t, double fp[], double fm[], double gp[],
                            double gm[],
                            fc2d_clawpack46_rpn2_t rpn2,
                            fc2d_clawpack46_rpt2_t rpt2,
                            fc2d_clawpack46_flux2_t flux2,
                            int block_corner_count[],int* ierror);

#define CLAWPACK46_FLUX2 FCLAW_F77_FUNC(clawpack46_flux2,CLAWPACK46_FLUX2)
void CLAWPACK46_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                      const int* maux,const int* mbc,const int* mx,
                      double q1d[], double dtdx1d[],
                      double aux1[], double aux2[], double aux3[],
                      double faddm[],double faddp[], double gaddm[],
                      double gaddp[],double cfl1d[], double fwave[],
                      double s[], double amdq[],double apdq[],double cqxx[],
                      double bmasdq[], double bpasdq[],
                      fc2d_clawpack46_rpn2_t rpn2,fc2d_clawpack46_rpt2_t rpt2,
                      const int* mwaves, const int* mcapa,
                      int method[], int mthlim[]);

#define CLAWPACK46_FLUX2FW FCLAW_F77_FUNC(clawpack46_flux2fw,CLAWPACK46_FLUX2FW)
void CLAWPACK46_FLUX2FW(const int* ixy,const int* maxm, const int* meqn, //
                        const int* maux,const int* mbc,const int* mx,
                        double q1d[], double dtdx1d[],
                        double aux1[], double aux2[], double aux3[],
                        double faddm[],double faddp[], double gaddm[],
                        double gaddp[],double cfl1d[], double fwave[],
                        double s[], double amdq[],double apdq[],double cqxx[],
                        double bmasdq[], double bpasdq[],
                        fc2d_clawpack46_rpn2_t rpn2,fc2d_clawpack46_rpt2_t rpt2,
                        const int* mwaves, const int* mcapa,
                        int method[], int mthlim[]);


#define CLAWPACK46_SET_BLOCK FCLAW_F77_FUNC(clawpack46_set_block,CLAWPACK46_SET_BLOCK)
void CLAWPACK46_SET_BLOCK(int* blockno);

#define FC2D_CLAWPACK46_GET_BLOCK FCLAW_F77_FUNC(fc2d_clawpack46_get_block, \
                                                 FC2D_CLAWPACK46_GET_BLOCK)
int FC2D_CLAWPACK46_GET_BLOCK();


#define CLAWPACK46_UNSET_BLOCK FCLAW_F77_FUNC(clawpack46_unset_block, \
                                              CLAWPACK46_UNSET_BLOCK)
void CLAWPACK46_UNSET_BLOCK();

/* -------------------------------------------------------------------------
   Old routines (that need to be replaced)
   ------------------------------------------------------------------------- */
void fc2d_clawpack46_define_auxarray(fclaw2d_domain_t* domain, ClawPatch *cp);
void fc2d_clawpack46_get_auxarray(fclaw2d_domain_t* domain,
                                 ClawPatch *cp, double **aux, int* maux);

void
    fc2d_clawpack46_set_options(fclaw2d_domain_t* domain,
                                fc2d_clawpack46_options_t* clawpack_options);

void
    fc2d_clawpack46_link_solvers(fclaw2d_domain_t* domain);

void
    fc2d_clawpack46_link_to_clawpatch();

/* -------------------------------------------------------------------------
   New routines
   ------------------------------------------------------------------------- */
void fc2d_clawpack46_define_auxarray2(fclaw2d_domain_t* domain, fclaw2d_patch_t* this_patch);

void fc2d_clawpack46_aux_data(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              double **aux, int* maux);

int fc2d_clawpack46_get_maux(fclaw2d_domain_t* domain);
void fc2d_clawpack46_maux(fclaw2d_domain_t* domain, int* maux);

void fc2d_clawpack46_register(fclaw_app_t* app, const char *configfile);

void fc2d_clawpack46_package_register(fclaw_app_t* app,
                                      fc2d_clawpack46_options_t* clawopt);

int fc2d_clawpack46_get_package_id();

fc2d_clawpack46_options_t* fc2d_clawpack46_get_options(fclaw2d_domain_t *domain);


/* -------------------------------------------------------------------------
   Routines that won't change
   ------------------------------------------------------------------------- */
void
    fc2d_clawpack46_setprob(fclaw2d_domain_t* domain);

void
    fc2d_clawpack46_setaux(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx);

void
    fc2d_clawpack46_qinit(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx);

void
    fc2d_clawpack46_b4step2(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt);

void
    fc2d_clawpack46_bc2(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        fclaw_bool intersects_bc[],
                        fclaw_bool time_interp);

void
    fc2d_clawpack46_src2(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt);


/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double
    fc2d_clawpack46_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

/* Use this ro return only the right hand side of the clawpack algorithm */
double
    fc2d_clawpack46_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);

double
fc2d_clawpack46_update(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif
