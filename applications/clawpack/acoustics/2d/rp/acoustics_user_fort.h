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

#ifndef ACOUSTICS_USER_H
#define ACOUSTICS_USER_H

#include <clawpack46_user_fort.h>
#include <clawpack5_user_fort.h>


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* ------------------------------ Clawpack 4.6 headers ------------------------------------------ */

#define CLAWPACK46_RPN2_MANIFOLD    FCLAW_F77_FUNC(clawpack46_rpn2_manifold,   \
                                                   CLAWPACK46_RPN2_MANIFOLD)
void CLAWPACK46_RPN2_MANIFOLD(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc,const int* mx, double ql[], double qr[],
                     double auxl[], double auxr[], double wave[],
                     double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2_MANIFOLD    FCLAW_F77_FUNC(clawpack46_rpt2_manifold,  \
                                                   CLAWPACK46_RPT2_MANIFOLD)
void CLAWPACK46_RPT2_MANIFOLD(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc, const int* mx, double ql[], double qr[],
                     double aux1[], double aux2[], double aux3[], const int* imp,
                     double dsdq[], double bmasdq[], double bpasdq[]);

#define CLAWPACK46_RPN2_ACOUSTICS_VC    FCLAW_F77_FUNC(clawpack46_rpn2_acoustics_vc,   \
                                                   CLAWPACK46_RPN2_ACOUSTICS_VC)
void CLAWPACK46_RPN2_ACOUSTICS_VC(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc,const int* mx, double ql[], double qr[],
                     double auxl[], double auxr[], double wave[],
                     double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2_ACOUSTICS_VC    FCLAW_F77_FUNC(clawpack46_rpt2_acoustics_vc,  \
                                                   CLAWPACK46_RPT2_ACOUSTICS_VC)
void CLAWPACK46_RPT2_ACOUSTICS_VC(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                     const int* mbc, const int* mx, double ql[], double qr[],
                     double aux1[], double aux2[], double aux3[], const int* imp,
                     double dsdq[], double bmasdq[], double bpasdq[]);

/* ------------------------------ Clawpack 4.6 headers ------------------------------------------ */

#define CLAWPACK5_RPN2_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpn2_manifold, \
                                               CLAWPACK5_RPN2_MANIFOLD)
void CLAWPACK5_RPN2_MANIFOLD(const int* ixy,const int* maxm, const int* meqn,
                                const int* mwaves, const int* maux,
                                const int* mbc,const int* mx,
                                double ql[], double qr[], double auxl[], double auxr[],
                                double wave[], double s[],double amdq[], double apdq[]);

#define CLAWPACK5_RPT2_MANIFOLD FCLAW_F77_FUNC(clawpack5_rpt2_manifold, \
                                               CLAWPACK5_RPT2_MANIFOLD)
void CLAWPACK5_RPT2_MANIFOLD(const int* ixy, const int* imp,
                             const int* maxm, const int* meqn,
                             const int* mwaves, const int* maux,
                             const int* mbc,const int* mx,
                             double ql[], double qr[],
                             double aux1[], double aux2[],
                             double aux3[],  double asdq[],
                             double bmasdq[], double bpasdq[]);

#define CLAWPACK5_RPN2_ACOUSTICS_VC FCLAW_F77_FUNC(clawpack5_rpn2_acoustics_vc, \
                                               CLAWPACK5_RPN2_ACOUSTICS_VC)
void CLAWPACK5_RPN2_ACOUSTICS_VC(const int* ixy,const int* maxm, const int* meqn,
                                const int* mwaves, const int* maux,
                                const int* mbc,const int* mx,
                                double ql[], double qr[], double auxl[], double auxr[],
                                double wave[], double s[],double amdq[], double apdq[]);

#define CLAWPACK5_RPT2_ACOUSTICS_VC FCLAW_F77_FUNC(clawpack5_rpt2_acoustics_vc, \
                                               CLAWPACK5_RPT2_ACOUSTICS_VC)
void CLAWPACK5_RPT2_ACOUSTICS_VC(const int* ixy, const int* imp,
                             const int* maxm, const int* meqn,
                             const int* mwaves, const int* maux,
                             const int* mbc,const int* mx,
                             double ql[], double qr[],
                             double aux1[], double aux2[],
                             double aux3[],  double asdq[],
                             double bmasdq[], double bpasdq[]);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
