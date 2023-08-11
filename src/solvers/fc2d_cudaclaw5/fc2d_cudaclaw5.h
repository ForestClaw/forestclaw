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

#ifndef FC2D_CUDACLAW5_H
#define FC2D_CUDACLAW5_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw_global;
struct fclaw2d_patch;

#if 0
    /* to get syntax highlighting */
#endif

/* --------------------------- Clawpack solver functions ------------------------------ */

/* Virtualize clawpack-specific wrapper functions */
typedef void (*cudaclaw5_src2_t)(struct fclaw_global* glob,
                                 struct fclaw2d_patch *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx,
                                 double t,
                                 double dt);

typedef void (*cudaclaw5_b4step2_t)(struct fclaw_global* glob,
                                    struct fclaw2d_patch *this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    double t,
                                    double dt);


/* -------------------------- Clawpack solver functions ------------------------------ */

typedef void (*cudaclaw5_fort_setprob_t)(void);

typedef void (*cudaclaw5_fort_bc2_t)(const int* meqn, const int* mbc,
                                     const int* mx, const int* my,
                                     const double* xlower, const double* ylower,
                                     const double* dx, const double* dy,
                                     const double q[], const int* maux,
                                     const double aux[], const double* t,
                                     const double* dt, const int mthbc[]);

typedef  void (*cudaclaw5_fort_qinit_t)(const int* meqn,const int* mbc,
                                        const int* mx, const int* my,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        double q[], const int* maux, double aux[]);

typedef void (*cudaclaw5_fort_setaux_t)(const int* mbc,
                                        const int* mx, const int* my,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        const int* maux, double aux[]);

typedef void (*cudaclaw5_fort_b4step2_t)(const int* mbc,
                                         const int* mx, const int* my, const int* meqn,
                                         double q[], const double* xlower,
                                         const double* ylower,
                                         const double* dx, const double* dy,
                                         const double* t, const double* dt,
                                         const int* maux, double aux[]);

typedef void (*cudaclaw5_fort_src2_t)(const int* meqn,
                                      const int* mbc, const int* mx,const int* my,
                                      const double* xlower, const double* ylower,
                                      const double* dx, const double* dy, double q[],
                                      const int* maux, double aux[], const double* t,
                                      const double* dt);


typedef void (*cudaclaw5_fort_rpn2_t)(const int* ixy,const int* maxm, const int* meqn,
                                      const int* mwaves, const int* maux,
                                      const int* mbc,const int* mx,
                                      double ql[], double qr[], double auxl[], double auxr[],
                                      double wave[], double s[],double amdq[], double apdq[]);

typedef void (*cudaclaw5_cuda_rpn2_t)(int idir, int meqn, int mwaves, int maux,
                                      double ql[], double qr[], 
                                      double auxl[], double auxr[],
                                      double wave[], double s[], 
                                      double amdq[], double apdq[]);

typedef void (*cudaclaw5_fort_rpt2_t)(const int* ixy, const int* imp, const int* maxm,
                                      const int* meqn,
                                       const int* mwaves, const int* maux, const int* mbc,
                                       const int* mx,
                                       double ql[], double qr[], double aux1[],
                                       double aux2[],
                                       double aux3[],  double asdq[],
                                       double bmasdq[], double bpasdq[]);


typedef void (*cudaclaw5_fort_flux2_t)(const int* ixy,const int* maxm, const int* meqn,
                                        const int* maux,const int* mbc,const int* mx,
                                        double q1d[], double dtdx1d[],
                                        double aux1[], double aux2[], double aux3[],
                                        double faddm[],double faddp[], double gaddm[],
                                        double gaddp[],double cfl1d[], double wave[],
                                        double s[], double amdq[],double apdq[],double cqxx[],
                                        double bmasdq[], double bpasdq[],
                                        cudaclaw5_fort_rpn2_t rpn2,
                                        cudaclaw5_fort_rpt2_t rpt2);


#if 0
typedef void (*cudaclaw5_fort_fluxfun_t)(const int* meqn, double q[], double aux[],
                                          double fq[]);
#endif

/* ------------------------------ CUDA functions ------------------------------------ */

double cudaclaw5_step2(struct fclaw_global* glob,
                       struct fclaw2d_patch* this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt);

/* ------------------------------------- Virtual table -------------------------------- */

typedef struct fc2d_cudaclaw5_vtable
{

    /* C wrappers */
    cudaclaw5_b4step2_t        b4step2;
    cudaclaw5_src2_t           src2;

    cudaclaw5_fort_setprob_t   fort_setprob;
    cudaclaw5_fort_qinit_t     fort_qinit;
    cudaclaw5_fort_setaux_t    fort_setaux;
    cudaclaw5_fort_bc2_t       fort_bc2;
    cudaclaw5_fort_b4step2_t   fort_b4step2;
    cudaclaw5_fort_src2_t      fort_src2;

    cudaclaw5_fort_rpn2_t      fort_rpn2;
    cudaclaw5_fort_rpt2_t      fort_rpt2;
    cudaclaw5_cuda_rpn2_t      cuda_rpn2;

    int is_set;
} fc2d_cudaclaw5_vtable_t;

void fc2d_cudaclaw5_solver_initialize();

fc2d_cudaclaw5_vtable_t* fc2d_cudaclaw5_vt();

/* ----------------------------- User access to solver functions ---------------------- */


void fc2d_cudaclaw5_setprob(struct fclaw_global* glob);

void fc2d_cudaclaw5_setaux(struct fclaw_global *glob,
                           struct fclaw2d_patch *this_patch,
                           int this_block_idx,
                           int this_patch_idx);

void fc2d_cudaclaw5_set_capacity(struct fclaw_global *glob,
                                 struct fclaw2d_patch *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx);

void fc2d_cudaclaw5_qinit(struct fclaw_global *glob,
                          struct fclaw2d_patch *this_patch,
                          int this_block_idx,
                          int this_patch_idx);

void fc2d_cudaclaw5_b4step2(struct fclaw_global *glob,
                            struct fclaw2d_patch *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt);

void fc2d_cudaclaw5_bc2(struct fclaw_global *glob,
                        struct fclaw2d_patch *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        int intersects_bc[],
                        int time_interp);

void fc2d_cudaclaw5_src2(struct fclaw_global *glob,
                         struct fclaw2d_patch *this_patch,
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


#endif /* !FC2D_CLAWPACH5_H */
