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

#ifndef TRANSPORT_USER_FORT_H
#define TRANSPORT_USER_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define USER46_B4STEP2_MANIFOLD FCLAW_F77_FUNC(user46_b4step2_manifold,USER46_B4STEP2_MANIFOLD)
void USER46_B4STEP2_MANIFOLD(const int* mx, const int* my, const int* mbc,
                             const double* dx, const double* dy,
                             const double* t, const int* maux, double aux[],
                             const int* blockno,
                             double xd[], double yd[], double zd[]);

#define USER5_B4STEP2_MANIFOLD FCLAW_F77_FUNC(user5_b4step2_manifold,USER5_B4STEP2_MANIFOLD)
void USER5_B4STEP2_MANIFOLD(const int* mx, const int* my, const int* mbc,
                            const double* dx, const double* dy,
                            const double* t, const int* maux, double aux[],
                            const int* blockno,
                            double xd[], double yd[], double zd[]);

#define USER_EXCEEDS_THRESHOLD FCLAW_F77_FUNC(user_exceeds_threshold, \
                                              USER_EXCEEDS_THRESHOLD)

int USER_EXCEEDS_THRESHOLD(int* blockno,
                           double qval[], 
                           double* qmin, double *qmax,
                           double quad[], 
                           double *dx, double *dy, 
                           double *xc, double *yc, 
                           int* tag_threshold, 
                           int* init_flag,
                           int* is_ghost);

#define CLAWPACK46_RPN2ADV FCLAW_F77_FUNC(clawpack46_rpn2adv,CLAWPACK46_RPN2ADV)
void CLAWPACK46_RPN2ADV(const int* ixy,const int* maxm, const int* meqn, 
                        const int* mwaves,
                        const int* mbc,const int* mx, 
                        double ql[], double qr[],
                        double auxl[], double auxr[], 
                        double wave[],
                        double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2ADV FCLAW_F77_FUNC(clawpack46_rpt2adv, CLAWPACK46_RPT2ADV)
void CLAWPACK46_RPT2ADV(const int* ixy, const int* maxm, const int* meqn, 
                        const int* mwaves,
                        const int* mbc, const int* mx, 
                        double ql[], double qr[],
                        double aux1[], double aux2[], 
                        double aux3[], const int* imp,
                        double dsdq[], double bmasdq[], double bpasdq[]);


#define CLAWPACK46_RPN2ADV_MANIFOLD FCLAW_F77_FUNC(clawpack46_rpn2adv_manifold, \
                                                   CLAWPACK46_RPN2ADV_MANIFOLD)
void CLAWPACK46_RPN2ADV_MANIFOLD(const int* ixy,const int* maxm, const int* meqn, 
                                 const int* mwaves,
                                 const int* mbc,const int* mx, 
                                 double ql[], double qr[],
                                 double auxl[], double auxr[], double wave[],
                                 double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2ADV_MANIFOLD FCLAW_F77_FUNC(clawpack46_rpt2adv_manifold, \
                                                   CLAWPACK46_RPT2ADV_MANIFOLD)
void CLAWPACK46_RPT2ADV_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, 
                                 const int* mwaves,
                                 const int* mbc, const int* mx, 
                                 double ql[], double qr[],
                                 double aux1[], double aux2[], double aux3[], 
                                 const int* imp,
                                 double dsdq[], double bmasdq[], double bpasdq[]);


#define CLAWPACK5_RPN2ADV FCLAW_F77_FUNC(clawpack5_rpn2adv,CLAWPACK5_RPN2ADV)
void CLAWPACK5_RPN2ADV(const int* ixy,const int* maxm, const int* meqn,
                       const int* mwaves, const int* maux,
                       const int* mbc,const int* mx,
                       double ql[], double qr[], double auxl[], double auxr[],
                       double wave[], double s[],double amdq[], double apdq[]);

#define CLAWPACK5_RPT2ADV FCLAW_F77_FUNC(clawpack5_rpt2adv, CLAWPACK5_RPT2ADV)
void CLAWPACK5_RPT2ADV(const int* ixy, const int* imp,
                       const int* maxm, const int* meqn,
                       const int* mwaves, const int* maux,
                       const int* mbc,const int* mx,
                       double ql[], double qr[],
                       double aux1[], double aux2[],
                       double aux3[],  double asdq[],
                       double bmasdq[], double bpasdq[]);

#define CLAWPACK5_RPN2ADV_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpn2adv_manifold,    \
                                                  CLAWPACK5_RPN2ADV_MANIFOLD)
void CLAWPACK5_RPN2ADV_MANIFOLD(const int* ixy,const int* maxm, const int* meqn,
                                const int* mwaves, const int* maux,
                                const int* mbc,const int* mx,
                                double ql[], double qr[], 
                                double auxl[], double auxr[],
                                double wave[], double s[],
                                double amdq[], double apdq[]);

#define CLAWPACK5_RPT2ADV_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpt2adv_manifold,    \
                                                  CLAWPACK5_RPT2ADV_MANIFOLD)
void CLAWPACK5_RPT2ADV_MANIFOLD(const int* ixy, const int* imp,
                                const int* maxm, const int* meqn,
                                const int* mwaves, const int* maux,
                                const int* mbc,const int* mx,
                                double ql[], double qr[],
                                double aux1[], double aux2[],
                                double aux3[],  double asdq[],
                                double bmasdq[], double bpasdq[]);

#define USER46_SETAUX_MANIFOLD FCLAW_F77_FUNC(user46_setaux_manifold, \
                                               USER46_SETAUX_MANIFOLD)

void USER46_SETAUX_MANIFOLD(const int* mbc,
                            const int* mx, const int* my,
                            const double* xlower, const double* ylower,
                            const double* dx, const double* dy,
                            const int* maux, double aux[],
                            const int* blockno,
                            double xd[], double yd[], double zd[],
                            double area[]);


#define USER5_SETAUX_MANIFOLD FCLAW_F77_FUNC(user5_setaux_manifold, \
                                             USER5_SETAUX_MANIFOLD)

void USER5_SETAUX_MANIFOLD(const int* mbc,
                           const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           const int* maux, double aux[],
                           const int* blockno,
                           double xd[], double yd[], double zd[],
                           double area[]);




#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
