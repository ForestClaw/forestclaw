/*
Copyright (c) 2019-2020 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef FC2D_THUNDEREGG_FORT_H
#define FC2D_THUNDEREGG_FORT_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */
#include "fc2d_thunderegg.h"    /* Needed for fc2d_thunderegg_fort_eval_bc_t def */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* - -------------------------------- BC functions ------------------------------------ */

#define THUNDEREGG_FORT_APPLY_BC_DEFAULT FCLAW_F77_FUNC(thunderegg_fort_apply_bc_default, \
                                                       THUNDEREGG_FORT_APPLY_BC_DEFAULT)

void THUNDEREGG_FORT_APPLY_BC_DEFAULT(const int* blockno, const  int* mx, const  int* my, 
                                     const  int* mbc, const  int* mfields, 
                                     const double* xlower, const double* ylower,
                                     const double* dx, const double* dy, 
                                     const double* t, 
                                     int intersects_bc[], int mthbc[], 
                                     double rhs[], fc2d_thunderegg_fort_eval_bc_t g_bc, 
                                     int *cons_check, double flux_sum[]);


#define THUNDEREGG_FORT_EVAL_BC_DEFAULT FCLAW_F77_FUNC(thunderegg_fort_eval_bc_default, \
                                                      THUNDEREGG_FORT_EVAL_BC_DEFAULT)

double THUNDEREGG_FORT_EVAL_BC_DEFAULT(const int* iface, const double* t, const double* x, const double* y);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif

