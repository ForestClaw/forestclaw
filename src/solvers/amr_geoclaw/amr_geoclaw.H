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

#ifndef AMR_AMR_GEOCLAW_H
#define AMR_AMR_GEOCLAW_H

#include "forestclaw2d.h"
#include "amr_options.h"
#include "clawpack_fort.H"

/* I know this is one 'extern C' too many, but I like the indentation here */
extern "C"
{
    void set_module_defaults();
    void setprob_();

    void bc2amr_(const double q[], const double aux[],
                 const int& nrol, const int& ncol, const int& meqn, const int& maux,
                 const int& hx, const int& hy, const int& level, const double& t,
                 const double& xlower, const double& ylower,
                 const double& xupper, const double& yupper,
                 const double& xperiodic, const double& yperiodic, const double spheredom);

    /* same as clawpack */
    void qinit_(const int& maxmx, const int& maxmy, const int& meqn,
                const int& mbc, const int& mx, const int& my,
                const double& xlower, const double& ylower,
                const double& dx, const double& dy,
                double q[], const int& maux, double aux[]);

    /* same as clawpack */
    void setaux_(const int& maxmx, const int& maxmy, const int& mbc,
                 const int& mx, const int& my,
                 const double& xlower, const double& ylower,
                 const double& dx, const double& dy,
                 const int& maux, double aux[]);

    /* Same as clawpack */
    void b4step2_(const int& maxmx, const int& maxmy, const int& mbc,
                  const int& mx, const int& my, const int& meqn,
                  double q[], const double& xlower, const double& ylower,
                  const double& dx, const double& dy,
                  const double& t, const double& dt,
                  const int& maux, double aux[]);

    /* same as clawpack */
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

typedef struct amr_geoclaw_parms amr_geoclaw_parms_t;
typedef struct amr_geoclaw_patch_data amr_geoclaw_patch_data_t;

/* Only one copy of geoclaw_parms for each run */
struct amr_geoclaw_parms
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

/* Reading in parameters */
void set_geoclaw_parms(fclaw2d_domain_t* domain,amr_geoclaw_parms_t* geoclaw_parms);
amr_geoclaw_parms_t* amr_geoclaw_parms_new(sc_options_t *opt);
void amr_geoclaw_checkparms(amr_geoclaw_parms_t *geoclaw_parms,amr_options_t* gparms);
void amr_geoclaw_parms_delete(amr_geoclaw_parms_t *geoclaw_parms);
void amr_geoclaw_postprocess_parms(amr_geoclaw_parms_t* geoclaw_parms);

/* Linking solvers */
void amr_geoclaw_link_solvers(fclaw2d_domain_t* domain);

/* Problem setup - user defined */
void amr_geoclaw_setprob(fclaw2d_domain_t* domain);

/* ----------------------------------------------------------------------
   Patch related parms
   ---------------------------------------------------------------------- */

/* Patch data is stored in each ClawPatch */
struct amr_geoclaw_patch_data
{
    FArrayBox auxarray;
    int maux;
};


void amr_geoclaw_link_to_clawpatch();


void amr_geoclaw_setaux(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx);

void amr_geoclaw_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx);

void amr_geoclaw_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

void amr_geoclaw_bc2(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      double t,
                      double dt,
                      fclaw_bool intersects_bc[]);


/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double amr_geoclaw_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

/* Use this ro return only the right hand side of the geoclaw algorithm */
double amr_geoclaw_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);


void dump_auxarray(fclaw2d_domain_t *domain, int dump_patchno);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif
