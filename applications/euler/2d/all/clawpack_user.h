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

#ifndef CLAWPACK_USER_H
#define CLAWPACK_USER_H

#include <fclaw2d_clawpatch.h>
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define CLAWPACK46_RPN2EU4 FCLAW_F77_FUNC(clawpack46_rpn2eu4,CLAWPACK46_RPN2EU4)
void CLAWPACK46_RPN2EU4(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                        const int* mbc,const int* mx, double ql[], double qr[],
                        double auxl[], double auxr[], double wave[],
                        double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2EU FCLAW_F77_FUNC(clawpack46_rpt2eu, CLAWPACK46_RPT2EU)
void CLAWPACK46_RPT2EU(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                        const int* mbc, const int* mx, double ql[], double qr[],
                        double aux1[], double aux2[], double aux3[], const int* imp,
                        double dsdq[], double bmasdq[], double bpasdq[]);

#define CLAWPACK5_RPN2EU4 FCLAW_F77_FUNC(clawpack5_rpn2eu4,CLAWPACK5_RPN2EU4)
void CLAWPACK5_RPN2EU4(const int* ixy,const int* maxm, const int* meqn,
                       const int* mwaves, const int* maux,
                       const int* mbc,const int* mx,
                       double ql[], double qr[], double auxl[], double auxr[],
                       double wave[], double s[],double amdq[], double apdq[]);


#define CLAWPACK5_RPT2EU FCLAW_F77_FUNC(clawpack5_rpt2eu, CLAWPACK5_RPT2EU)
void CLAWPACK5_RPT2EU(const int* ixy, const int* imp,
                       const int* maxm, const int* meqn,
                       const int* mwaves, const int* maux,
                       const int* mbc,const int* mx,
                       double ql[], double qr[],
                       double aux1[], double aux2[],
                       double aux3[],  double asdq[],
                       double bmasdq[], double bpasdq[]);

#define CLAWPACK5_RPN2EU4_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpn2eu4_manifold,    \
                                                  CLAWPACK5_RPN2EU4_MANIFOLD)
void CLAWPACK5_RPN2EU4_MANIFOLD(const int* ixy,const int* maxm, const int* meqn,
                                const int* mwaves, const int* maux,
                                const int* mbc,const int* mx,
                                double ql[], double qr[], double auxl[], double auxr[],
                                double wave[], double s[],double amdq[], double apdq[]);


#define CLAWPACK5_RPT2EU_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpt2eu_manifold, \
                                                  CLAWPACK5_RPT2EU_MANIFOLD)
void CLAWPACK5_RPT2EU_MANIFOLD(const int* ixy, const int* imp,
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
