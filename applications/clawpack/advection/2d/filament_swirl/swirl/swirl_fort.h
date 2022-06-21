
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

#ifndef SWIRL_FORT_H
#define SWIRL_FORT_H

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


#define SWIRL_SETPROB            FCLAW_F77_FUNC(swirl_setprob,           SWIRL_SETPROB)
void SWIRL_SETPROB();

#define SWIRL_CLAWPACK46_QINIT   FCLAW_F77_FUNC(swirl_clawpack46_qinit,  SWIRL_CLAWPACK46_QINIT)
void SWIRL_CLAWPACK46_QINIT(const int* maxmx, const int* maxmy, const int* meqn,
                      const int* mbc, const int* mx, const int* my,
                      const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      double q[], const int* maux, double aux[]);

#define SWIRL_CLAWPACK5_QINIT   FCLAW_F77_FUNC(swirl_clawpack5_qinit,   SWIRL_CLAWPACK5_QINIT)
void SWIRL_CLAWPACK5_QINIT(const int* meqn,const int* mbc,
                     const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy,
                     double q[], const int* maux, double aux[]);

#define SWIRL_CLAWPACK46_B4STEP2 FCLAW_F77_FUNC(swirl_clawpack46_b4step2,SWIRL_CLAWPACK46_B4STEP2)
void SWIRL_CLAWPACK46_B4STEP2(const int* maxmx, const int* maxmy, const int* mbc,
                        const int* mx, const int* my, const int* meqn,
                        double q[], const double* xlower, const double* ylower,
                        const double* dx, const double* dy,
                        const double* t, const double* dt,
                        const int* maux, double aux[]);

#define SWIRL_CLAWPACK5_B4STEP2 FCLAW_F77_FUNC(swirl_clawpack5_b4step2, SWIRL_CLAWPACK5_B4STEP2)
void SWIRL_CLAWPACK5_B4STEP2(const int* mbc,
                       const int* mx, const int* my, const int* meqn,
                       double q[], const double* xlower,
                       const double* ylower,
                       const double* dx, const double* dy,
                       const double* t, const double* dt,
                       const int* maux, double aux[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
