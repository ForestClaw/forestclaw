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

typedef void (*fc2d_clawpack46_setprob_t)();

typedef void (*fc2d_clawpack46_bc2_t)(const int& maxmx, const int& maxmy,
                                      const int& meqn, const int& mbc,
                                      const int& mx, const int& my,
                                      const double& xlower, const double& ylower,
                                      const double& dx, const double& dy,
                                      const double q[], const int& maux,
                                      const double aux[], const double& t,
                                      const double& dt, const int mthbc[]);

typedef  void (*fc2d_clawpack46_qinit_t)(const int& maxmx, const int& maxmy,
                                         const int& meqn,const int& mbc,
                                         const int& mx, const int& my,
                                         const double& xlower, const double& ylower,
                                         const double& dx, const double& dy,
                                         double q[], const int& maux, double aux[]);

typedef void (*fc2d_clawpack46_setaux_t)(const int& maxmx, const int& maxmy, const int& mbc,
                                         const int& mx, const int& my,
                                         const double& xlower, const double& ylower,
                                         const double& dx, const double& dy,
                                         const int& maux, double aux[]);

typedef void (*fc2d_clawpack46_b4step2_t)(const int& maxmx, const int& maxmy,
                                          const int& mbc,
                                          const int& mx, const int& my, const int& meqn,
                                          double q[], const double& xlower,
                                          const double& ylower,
                                          const double& dx, const double& dy,
                                          const double& t, const double& dt,
                                          const int& maux, double aux[]);

typedef void (*fc2d_clawpack46_src2_t)(const int& maxmx, const int& maxmy, const int& meqn,
                                       const int& mbc, const int& mx,const int& my,
                                       const double& xlower, const double& ylower,
                                       const double& dx, const double& dy, double q[],
                                       const int& maux, double aux[], const double& t,
                                       const double& dt);

typedef struct fc2d_clawpack46_vtable
{
    fc2d_clawpack46_setprob_t setprob;
    fc2d_clawpack46_bc2_t bc2;
    fc2d_clawpack46_qinit_t qinit;
    fc2d_clawpack46_setaux_t setaux;
    fc2d_clawpack46_b4step2_t b4step2;
    fc2d_clawpack46_src2_t src2;
} fc2d_clawpack46_vtable_t;

void fc2d_clawpack46_set_vtable(const fc2d_clawpack46_vtable_t* vt);


extern "C"
{
    void clawpack46_bc2_(const int& maxmx, const int& maxmy, const int& meqn,
                         const int& mbc, const int& mx, const int& my,
                         const double& xlower, const double& ylower,
                         const double& dx, const double& dy, const double q[],
                         const int& maux, const double aux[], const double& t,
                         const double& dt, const int mthbc[]);

    /* --------------------------------------------------------------------
       Classic routines
       - These are provided only for convenience;  these files are not
         compiled into the library.
       -------------------------------------------------------------------- */
    void setprob_();

    void qinit_(const int& maxmx, const int& maxmy, const int& meqn,
                const int& mbc, const int& mx, const int& my,
                const double& xlower, const double& ylower,
                const double& dx, const double& dy,
                double q[], const int& maux, double aux[]);

    void setaux_(const int& maxmx, const int& maxmy, const int& mbc,
                 const int& mx, const int& my,
                 const double& xlower, const double& ylower,
                 const double& dx, const double& dy,
                 const int& maux, double aux[]);

    void b4step2_(const int& maxmx, const int& maxmy, const int& mbc,
                  const int& mx, const int& my, const int& meqn,
                  double q[], const double& xlower, const double& ylower,
                  const double& dx, const double& dy,
                  const double& t, const double& dt,
                  const int& maux, double aux[]);

    void src2_(const int& maxmx, const int& maxmy, const int& meqn,
               const int& mbc, const int& mx,const int& my,
               const double& xlower, const double& ylower,
               const double& dx, const double& dy, double q[],
               const int& maux, double aux[], const double& t,
               const double& dt);

    /* --------------------------------------------------------------------
       Time stepping
       -------------------------------------------------------------------- */
    void clawpatch2_(const int& maxm, const int& meqn, const int& maux, const int& mbc,
                     const int method[], const int mthlim[], const int& mcapa,
                     const int& mwaves, const int& mx, const int& my, double qold[],
                     double auxold[], const double& dx, const double& dy, const double& dt,
                     const double& cfl, double work[], const int& mwork, const double& xlower,
                     const double& ylower, const int& level, const double& t, double fp[],
                     double fm[], double gp[], double gm[]);

}
#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define SET_CORNERS FCLAW_F77_FUNC(set_corners, SET_CORNERS)
void SET_CORNERS(int block_corners_count[]);

#define CLAWPACK_SET_BLOCK FCLAW_F77_FUNC(clawpack_set_block, CLAWPACK_SET_BLOCK)
void CLAWPACK_SET_BLOCK(int* blockno);

#define CLAWPACK_GET_BLOCK FCLAW_F77_FUNC(clawpack_get_block, CLAWPACK_GET_BLOCK)
int CLAWPACK_GET_BLOCK();


void
fc2d_clawpack46_define_auxarray(fclaw2d_domain_t* domain,
                                ClawPatch *cp);

void
fc2d_clawpack46_get_auxarray(fclaw2d_domain_t* domain,
                             ClawPatch *cp, double **aux, int* maux);

void
fc2d_clawpack46_set_options(fclaw2d_domain_t* domain,
                            fc2d_clawpack46_options_t* clawpack_options);

void
fc2d_clawpack46_link_solvers(fclaw2d_domain_t* domain);

void
fc2d_clawpack46_setprob();

void
fc2d_clawpack46_link_to_clawpatch();


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
