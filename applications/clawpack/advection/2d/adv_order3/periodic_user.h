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

#ifndef PERIODIC_USER_H
#define PERIODIC_USER_H

#include <fclaw2d_include_all.h>

#include <fc2d_clawpack46.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int claw_version;
    int is_registered;

} user_options_t;

void periodic_link_solvers(fclaw2d_global_t *glob);

void periodic_problem_setup(fclaw2d_global_t* glob);

/* ------------------------------------- Options ---------------------------------------*/
user_options_t* periodic_options_register (fclaw_app_t * app,
                                        const char *configfile);

void periodic_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* periodic_get_options(fclaw2d_global_t* glob);

/* ------------------------------------ Fortran ----------------------------------------*/
#define PERIODIC_SETPROB FCLAW_F77_FUNC(periodic_setprob, PERIODIC_SETPROB)
void PERIODIC_SETPROB();

#define PERIODIC_FLUX2 FCLAW_F77_FUNC(periodic_flux2,PERIODIC_FLUX2)
void PERIODIC_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                      const int* maux,const int* mbc,const int* mx,
                      double q1d[], double dtdx1d[],
                      double aux1[], double aux2[], double aux3[],
                      double faddm[],double faddp[], double gaddm[],
                      double gaddp[],double cfl1d[], double fwave[],
                      double s[], double amdq[],double apdq[],double cqxx[],
                      double bmasdq[], double bpasdq[],
                      clawpack46_fort_rpn2_t rpn2, clawpack46_fort_rpt2_t rpt2,
                      const int* mwaves, const int* mcapa,
                      int method[], int mthlim[]);


#define PERIODIC_FORT_INTERPOLATE2FINE FCLAW_F77_FUNC( \
            periodic_fort_interpolate2fine, PERIODIC_FORT_INTERPOLATE2FINE)
void PERIODIC_FORT_INTERPOLATE2FINE(const int* mx,const int* my,const int* mbc,
                                    const int* meqn, double qcoarse[], double qfine[],
                                    double areacoarse[], double areafine[],
                                    const int* igrid,const int* manifold);

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





#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
