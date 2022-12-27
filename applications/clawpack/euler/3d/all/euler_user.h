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

#ifndef EULER_USER_3D_H
#define EULER_USER_3D_H

#include <fclaw2d_include_all.h>


/* Headers for both Clawpack 4.6 and  Clawpack 5.0 */
#include <fclaw3dx_clawpatch.h>
#include <fclaw3dx_clawpatch_options.h>
#include <fclaw3dx_clawpatch_fort.h>


/* Clawpack 4.6 headers */  
#include <fc3d_clawpack46.h>  
#include <fc3d_clawpack46_options.h>
#include <fc3d_clawpack46_fort.h>  
#include <fc3d_clawpack46_user_fort.h>  
#include <fclaw3dx_clawpatch46_fort.h>

#if 0
/* Clawpack 5.0 headers */
#include <fc2d_clawpack5.h>
#include <fc2d_clawpack5_options.h>
#include <fc2d_clawpack5_fort.h>
#include <clawpack5_user_fort.h>
#include <fclaw3dx_clawpatch5_fort.h>
#endif


/* Headers for common FORTRAN files */
//#include "advection_user_fort3.h"

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax */
#endif


#define EULER3D_SETAUX_MANIFOLD FCLAW_F77_FUNC(euler3d_setaux_manifold, \
                                               EULER3D_SETAUX_MANIFOLD)

void EULER3D_SETAUX_MANIFOLD(const int* mbc,
                            const int* mx, const int* my, const int*mz,
                            const int* mcapa, 
                            const double* xlower, const double* ylower,
                            const double* zlower,
                            const double* dx, const double* dy,
                            const double *dz,
                            const int* maux, double aux[],
                            const int* blockno,
                            double xrot[], double yrot[], double zrot[],
                            double volume[],double faceareas[]);


#define CLAWPACK46_RPN3_MAPPED FCLAW_F77_FUNC(clawpack46_rpn3_mapped,  CLAWPACK46_RPN3_MAPPED)
void CLAWPACK46_RPN3_MAPPED(const int* ixyz,const int* maxm, 
                     const int* meqn, const int* mwaves,
                     const int* maux, const int* mbc, const int* mx, 
                     double ql[], double qr[],
                     double auxl[], double auxr[], double wave[],
                     double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT3_MAPPED FCLAW_F77_FUNC(clawpack46_rpt3_mapped, CLAWPACK46_RPT3_MAPPED)
void CLAWPACK46_RPT3_MAPPED(const int* ixyz, const int* icoor, const int* imp,
                     const int *maxm, const int* meqn, const int* mwaves, 
                     const int *maux, 
                     const int* mbc, const int* mx, double ql[], double qr[],
                     double aux1[], double aux2[], double aux3[], 
                     double asdq[], double bmasdq[], double bpasdq[]);

#define CLAWPACK46_RPTT3_MAPPED    FCLAW_F77_FUNC(clawpack46_rptt3_mapped, CLAWPACK46_RPTT3_MAPPED)
void CLAWPACK46_RPTT3_MAPPED(const int* ixyz, const int* icoor, const int* imp,
                      const int* impt, const int* maxm, const int* meqn,
                      const int* mwaves, const int* maux,
                      const int* mbc,const int* mx,
                      double ql[], double qr[],
                      double aux1[], double aux2[],
                      double aux3[],  double bsasdq[],
                      double cmbsasdq[], double cpbsasdq[]);



#ifdef __cplusplus
}
#endif

#endif
