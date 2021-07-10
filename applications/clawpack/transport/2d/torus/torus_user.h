/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef TORUS_USER_H
#define TORUS_USER_H

#include "../all/transport_user.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
#endif

/* --------------------------
   Headers for both versions
   -------------------------- */

typedef struct user_options
{
    int example;    

    double alpha;     /* Ratio of inner radius to outer radius */
    double beta;
    double revs_per_s;

    int claw_version;

    int is_registered;
}
user_options_t;

void torus_link_solvers(fclaw2d_global_t *glob);

user_options_t* torus_options_register (fclaw_app_t * app,
                                       const char *configfile);

void torus_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* torus_get_options(fclaw2d_global_t* glob);


fclaw2d_map_context_t *
    fclaw2d_map_new_torus (fclaw2d_map_context_t* brick,
                           const double scale[],
                           const double rotate[],
                           const double alpha,
                           const double beta);

/* ----------------------
   Clawpack 4.6 headers
   ---------------------- */
#define TORUS46_SETAUX  FCLAW_F77_FUNC(torus46_setaux, TORUS46_SETAUX)
void TORUS46_SETAUX(const int* mbc, const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* maux, double aux[], const int* blockno,
                    double area[], double edgelengths[], 
                    double xnormals[], double ynormals[], 
                    double surfnormals[]);

#define RPN2CONS_FW_MANIFOLD FCLAW_F77_FUNC(rpn2cons_fw_manifold, RPN2CONS_FW_MANIFOLD)
void RPN2CONS_FW_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                          const int* mwaves,
                          const int* mbc, const int* mx, double ql[], double qr[],
                          double auxl[], double auxr[], double fwave[],
                          double s[], double amdq[], double apdq[], const int* maux);


#define RPT2CONS_MANIFOLD FCLAW_F77_FUNC(rpt2cons_manifold, RPT2CONS_MANIFOLD)
void RPT2CONS_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                       const int* mwaves, const int* maux,
                       const int* mbc, const int* mx, double ql[], double qr[],
                       double aux1[], double aux2[], double aux3[], const int* imp,
                       double dsdq[], double bmasdq[], double bpasdq[]);


#define RPN2_CONS_UPDATE_MANIFOLD FCLAW_F77_FUNC(rpn2_cons_update_manifold, \
                                                 RPN2_CONS_UPDATE_MANIFOLD)

void RPN2_CONS_UPDATE_MANIFOLD(const int* meqn, const int* maux, const int* idir,
                               const int* iface, double q[], 
                               double aux_center[], double aux_edge[],
                               double flux[]);

#define RPN2_CONS_UPDATE_ZERO FCLAW_F77_FUNC(rpn2_cons_update_zero, \
                                             RPN2_CONS_UPDATE_ZERO)

void RPN2_CONS_UPDATE_ZERO(const int* meqn, const int* maux, const int* idir,
                           const int* iface,
                           double q[], double aux_center[], double aux_edge[],
                           double flux[]);


/* ----------------------
   Clawpack 5.x headers
   ---------------------- */
#if 0
#define TORUS5_SETAUX  FCLAW_F77_FUNC(torus5_setaux,  TORUS5_SETAUX)
void TORUS5_SETAUX(const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   const int* maux, double aux[]);
#endif                   



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !TORUS_USER_H */
