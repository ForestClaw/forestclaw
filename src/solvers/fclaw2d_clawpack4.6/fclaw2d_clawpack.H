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

#ifndef FCLAW2D_CLAWPACK_H
#define FCLAW2D_CLAWPACK_H

#include "forestclaw2d.h"
#include "amr_options.h"
#include "clawpack_fort.H"
#include "amr_forestclaw.H"
#include "fclaw2d_options.h"


extern "C"
{
    void setprob_();

    void bc2_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
              const int& mx, const int& my, const double& xlower, const double& ylower,
              const double& dx, const double& dy, const double q[], const int& maux,
              const double aux[], const double& t, const double& dt, const int mthbc[]);

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

typedef struct fclaw2d_clawpack_parms fclaw2d_clawpack_parms_t;
typedef struct fclaw2d_clawpack_patch_data fclaw2d_clawpack_patch_data_t;

/* Only one copy of fclaw2d_clawpack_parms for each run */
struct fclaw2d_clawpack_parms
{
    int mwaves;
    int maux;

    const char *order_string;
    int *order;   /* order[2] */

    const char *mthlim_string;
    int *mthlim;  /* mthlim[mwaves] */

    int method[7];  /* method[7] */
    int mcapa;
    int src_term;
    void *user;

};


/* ----------------------------------------------------------------------
   Routines that get called once when setting up the problem and solvers.
   ---------------------------------------------------------------------- */
/* clawpack specific arguments */
void fclaw2d_clawpack_define_auxarray(fclaw2d_domain_t* domain, ClawPatch *cp);
void fclaw2d_clawpack_get_auxarray(fclaw2d_domain_t* domain,
                               ClawPatch *cp, double **aux, int* maux);

/* Reading in parameters */
void set_clawpack_parms(fclaw2d_domain_t* domain,fclaw2d_clawpack_parms_t* clawpack_parms);
fclaw2d_clawpack_parms_t* fclaw2d_clawpack_parms_new(sc_options_t *opt);
void fclaw2d_clawpack_checkparms(fclaw2d_clawpack_parms_t *clawpack_parms,
                                 amr_options_t* gparms);

void fclaw2d_clawpack_parms_delete(fclaw2d_clawpack_parms_t *clawpack_parms);
void fclaw2d_clawpack_postprocess_parms(fclaw2d_clawpack_parms_t* clawpack_parms);

fclaw2d_clawpack_parms_t*  clawpack46_new_options();
void clawpack46_read_options_from_file(sc_options_t* opt);
int clawpack46_options_read_from_file(sc_options_t* opt,int log_priority);

void clawpack46_register_options (sc_options_t* opt, fclaw2d_clawpack_parms_t* clawpack_parms);
void clawpack46_postprocess_parms(fclaw2d_clawpack_parms_t* clawpack_parms);
int clawpack46_checkparms(sc_options_t *opt,
                           fclaw2d_clawpack_parms_t *clawpack_parms,
                           amr_options_t* gparms,
                           int lp);


/* Linking solvers */
void fclaw2d_clawpack_link_solvers(fclaw2d_domain_t* domain);

/* Problem setup - user defined */
void fclaw2d_clawpack_setprob(fclaw2d_domain_t* domain);

/* ----------------------------------------------------------------------
   Patch related parms
   ---------------------------------------------------------------------- */

/* Patch data is stored in each ClawPatch */
struct fclaw2d_clawpack_patch_data
{
    FArrayBox auxarray;
    int maux;
};


void fclaw2d_clawpack_link_to_clawpatch();


void fclaw2d_clawpack_setaux(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx);

void fclaw2d_clawpack_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx);

void fclaw2d_clawpack_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

void fclaw2d_clawpack_bc2(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      double t,
                      double dt,
                      fclaw_bool intersects_bc[],
                      fclaw_bool time_interp);

void fclaw2d_clawpack_src2(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt);


/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double fclaw2d_clawpack_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

/* Use this ro return only the right hand side of the clawpack algorithm */
double fclaw2d_clawpack_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);

double fclaw2d_clawpack_update(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt);

void dump_auxarray(fclaw2d_domain_t *domain, int dump_patchno);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif
