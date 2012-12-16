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

#ifndef AMR_AMR_WAVEPROP_H
#define AMR_AMR_WAVEPROP_H

#include "forestclaw2d.h"
#include "amr_options.h"

/* I know this is one 'extern C' too many, but I like the indentation here */
extern "C"
{
    void bc2_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
              const int& mx, const int& my, const double& xlower, const double& ylower,
              const double& dx, const double& dy, const double q[], const int& maux,
              const double aux[], const double& t, const double& dt, const int mthbc[]);

    void bc2_manifold_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
                       const int& mx, const int& my, const double& xlower, const double& ylower,
                       const double& dx, const double& dy, const double q[], const int& maux,
                       const double aux[], const double& t, const double& dt, const int mthbc[],
                       double xp[], double yp[], double zp[],
                       double xd[],double yd[],double zd[],
                       const double xface_normals[], const double yface_normals[]);

    void qinit_(const int& maxmx, const int& maxmy, const int& meqn,
                const int& mbc, const int& mx, const int& my,
                const double& xlower, const double& ylower,
                const double& dx, const double& dy,
                double q[], const int& maux, double aux[]);

    void qinit_manifold_(const int& maxmx, const int& maxmy, const int& meqn,
                         const int& mbc, const int& mx, const int& my,
                         const double& xlower, const double& ylower,
                         const double& dx, const double& dy,
                         double q[], const int& maux, double aux[],
                         double xp[], double yp[], double zp[],
                         double xd[], double yd[], double zd[],
                         const int& blockno);

    void setaux_(const int& maxmx, const int& maxmy, const int& mbc,
                 const int& mx, const int& my,
                 const double& xlower, const double& ylower,
                 const double& dx, const double& dy,
                 const int& maux, double aux[]);

    void setaux_manifold_(const int& maxmx, const int& maxmy, const int& mbc,
                          const int& mx, const int& my,
                          const double& xlower, const double& ylower,
                          const double& dx, const double& dy,
                          const int& maux, double aux[],
                          double xp[], double yp[], double zp[],
                          double xd[], double yd[], double zd[],double area[]);

    void b4step2_(const int& maxmx, const int& maxmy, const int& mbc,
                  const int& mx, const int& my, const int& meqn,
                  double q[], const double& xlower, const double& ylower,
                  const double& dx, const double& dy,
                  const double& t, const double& dt,
                  const int& maux, double aux[]);


    void b4step2_manifold_(const int& maxmx, const int& maxmy, const int& mbc,
                          const int& mx, const int& my, const int& meqn,
                          double q[], const double& xlower, const double& ylower,
                          const double& dx, const double& dy,
                          const double& t, const double& dt,
                          const int& maux, double aux[],
                          double xp[], double yp[], double zp[],
                          double xd[], double yd[], double zd[]);

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

typedef struct amr_waveprop_parms amr_waveprop_parms_t;

struct amr_waveprop_parms
{
    int mwaves;
    int maux;

    const char *order_string;
    int *order;   /* order[2] */

    const char *mthlim_string;
    int *mthlim;  /* mthlim[mwaves] */

    int method[7];  /* method[7] */
    int mcapa;
    double max_cfl;
    double desired_cfl;
    int src_term;
    void *user;

};

typedef struct amr_waveprop_patch_data
{
    FArrayBox auxarray;
    int maux;
} amr_waveprop_patch_data_t;

void amr_waveprop_readparms(sc_options_t *opt,  amr_options_t *gparms);
void amr_waveprop_checkparms(amr_options_t *gparms);

void amr_waveprop_setaux(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx);

void amr_waveprop_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx);

void amr_waveprop_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);


/* Use this ro return only the right hand side of the waveprop algorithm */
double amr_waveprop_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);

/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double amr_waveprop_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

void amr_waveprop_bc2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt,
                          fclaw_bool intersects_bc[]);


void dump_auxarray(fclaw2d_domain_t *domain, int dump_patchno);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif
