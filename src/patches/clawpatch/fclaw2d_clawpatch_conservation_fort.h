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

#ifndef FCLAW2D_CLAWPATCH_CONSERVATION_FORT_H
#define FCLAW2D_CLAWPATCH_CONSERVATION_FORT_H

#include <fclaw_base.h>

/**
 * @file
 * C declarations of Fortran subroutines
 */
#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
/* To fix syntax highlighting below */
#endif    

/** Fortran subroutine name */
#define CLAWPATCH_TIME_SYNC_SETUP FCLAW_F77_FUNC(clawpatch_time_sync_setup, \
                                                 CLAWPATCH_TIME_SYNC_SETUP)

/** @copydoc clawpatch_time_sync_setup() */
void CLAWPATCH_TIME_SYNC_SETUP(const int* mx,const int *my, const int* mbc, 
                                  const double* dx, const double *dy,
                                  double area[], double edge_lengths[],
                                  double area0[], double area1[], 
                                  double area2[], double area3[],
                                  double el0[], double el1[], 
                                  double el2[], double el3[], 
                                  const int* manifold);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif

