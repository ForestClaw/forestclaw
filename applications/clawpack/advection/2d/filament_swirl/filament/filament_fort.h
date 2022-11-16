
/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FILAMENT_FORT_H
#define FILAMENT_FORT_H

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


/* --------------------------------------------------------------------------------
   Clawpack 5.0 routines

   These are provided for user convenience.  These files are not compiled
   into the library, but should be provided by the user.

   These signatures can be used if the user file matches these signatures 
   and subroutine name. Otherwise, the user should provide their own headers.
   ------------------------------------------------------------------------------- */

/* Macros for C/Fortran portability */

#define FILAMENT_SETPROB            FCLAW_F77_FUNC(filament_setprob,           FILAMENT_SETPROB)
void FILAMENT_SETPROB();

#define FILAMENT_CLAWPACK46_QINIT   FCLAW_F77_FUNC(filament_clawpack46_qinit,  FILAMENT_CLAWPACK46_QINIT)
void FILAMENT_CLAWPACK46_QINIT(const int* maxmx, const int* maxmy, const int* meqn,
                      const int* mbc, const int* mx, const int* my,
                      const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      double q[], const int* maux, double aux[]);

#define FILAMENT_CLAWPACK5_QINIT   FCLAW_F77_FUNC(filament_clawpack5_qinit,   FILAMENT_CLAWPACK5_QINIT)
void FILAMENT_CLAWPACK5_QINIT(const int* meqn,const int* mbc,
                     const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy,
                     double q[], const int* maux, double aux[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
