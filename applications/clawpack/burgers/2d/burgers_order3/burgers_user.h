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

#ifndef BURGERS_USER_H
#define BURGERS_USER_H

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
    int example; 
    
    int claw_version;
    int is_registered;

} user_options_t;

void burgers_link_solvers(fclaw_global_t *glob);

/* ------------------------------------- Options ---------------------------------------*/
user_options_t* burgers_options_register (fclaw_app_t * app,
                                        const char *configfile);

void burgers_options_store (fclaw_global_t* glob, user_options_t* user);

const user_options_t* burgers_get_options(fclaw_global_t* glob);

/* ------------------------------------ Fortran ----------------------------------------*/
#define BURGERS_FLUX2 FCLAW_F77_FUNC(burgers_flux2,BURGERS_FLUX2)
void BURGERS_FLUX2(const int* ixy,const int* maxm, const int* meqn,
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

#define RPN2CONS_UPDATE FCLAW_F77_FUNC(rpn2cons_update,RPN2CONS_UPDATE)

void RPN2CONS_UPDATE(const int* meqn, const int* maux, const int* idir, const int* iface,
                      double q[], double aux_center[], double aux_edge[], double flux[]);


#define RPN2CONS_UPDATE_ORDER2 FCLAW_F77_FUNC(rpn2cons_update_order2, \
                                              RPN2CONS_UPDATE_ORDER2)

void RPN2CONS_UPDATE_ORDER2(const int* meqn, const int* maux, 
                            const int* idir, const int* iface,
                            double q[], double aux_center[], 
                            double aux_edge[], double flux[]);


#define CLAWPACK46_RPN2BU    FCLAW_F77_FUNC(clawpack46_rpn2bu,CLAWPACK46_RPN2BU)
void CLAWPACK46_RPN2BU(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc,const int* mx, double ql[], double qr[],
                     double auxl[], double auxr[], double wave[],
                     double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2BU    FCLAW_F77_FUNC(clawpack46_rpt2bu, CLAWPACK46_RPT2BU)
void CLAWPACK46_RPT2BU(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                       const int* mbc, const int* mx, double ql[], double qr[],
                       double aux1[], double aux2[], double aux3[], const int* imp,
                       double dsdq[], double bmasdq[], double bpasdq[]);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
