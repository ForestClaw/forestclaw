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

#ifndef FCLAW2D_CLAWPATCH_FORT_H
#define FCLAW2D_CLAWPATCH_FORT_H

#include <fclaw2d_defs.h>

#if FCLAW2D_PATCHDIM == 2
#include "fclaw2d_clawpatch_fort2.h"
#else
#include "fclaw2d_clawpatch_fort3.h"
#endif


#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

struct fclaw2d_patch_transform_data;  /* Should be replaced by long int?  */

/* Functions defined here are implemented in individual solvers (clawpack 4.6 and 
   clawpack 5.0) */

/* --------------------- Dimension independent defs and routines ---------------------- */

/* These headers are independent of dimension and clawpack version */
typedef void  (*clawpatch_fort_header_ascii_t)(const char* matname1,const char* matname2,
                                               const double* time, const int* meqn, 
                                               const int* maux, const int* ngrids);

/* Even though this is for 3d patches, we still assume that tagging can only be
   dependent on two dimensions */
typedef int (*clawpatch_fort_exceeds_threshold_t)(const int *blockno,
                                                  const double *qval, 
                                                  const double *qmin, 
                                                  const double *qmax,
                                                  const double quad[], 
                                                  const double *dx, 
                                                  const double *dy, 
                                                  const double *xc, 
                                                  const double *yc, 
                                                  const double *tag_threshold,
                                                  const int    *init_flag,
                                                  const int    *is_ghost);

/* ----------------------------- Fortran headers ---------------------------------------*/

#define FCLAW2D_CLAWPATCH_GET_REFINEMENT_CRITERIA \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_get_refinement_criteria, \
                                 FCLAW2D_CLAWPATCH_GET_REFINEMENT_CRITERIA)
int FCLAW2D_CLAWPATCH_GET_REFINEMENT_CRITERIA();


/* ------------------------------- General threshold ---------------------------------- */
#define FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_exceeds_threshold, \
                                  FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD)

int FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD(const int *blockno,
                                        const double *qval, 
                                        const double *qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *xc, 
                                        const double *yc, 
                                        const double *tag_threshold,
                                        const int *init_flag,
                                        const int *is_ghost);

/* ----------------------------- Value threshold -------------------------------------- */
#define FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_value_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH)
    
int FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH(const int* blockno,
                                       const double *qval, 
                                       const double* qmin, 
                                       const double *qmax,
                                       const double quad[], 
                                       const double *dx, 
                                       const double *dy, 
                                       const double *xc, 
                                       const double *yc, 
                                       const double* tag_threshold,
                                       const int* init_flag,
                                       const int* is_ghost);
    
/* ----------------------------- difference threshold --------------------------------- */

#define FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_difference_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH)

int FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH(const int    *blockno,
                                            const double *qval, 
                                            const double *qmin, 
                                            const double *qmax,
                                            const double quad[], 
                                            const double *dx, 
                                            const double *dy, 
                                            const double *xc, 
                                            const double *yc, 
                                            const double *tag_threshold,
                                            const int *init_flag,
                                            const int *is_ghost);

/* --------------------------------- minmax threshold --------------------------------- */

#define FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_minmax_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH)

int FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH(const int *blockno,
                                        const double *qval, 
                                        const double* qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *xc, 
                                        const double *yc, 
                                        const double *tag_threshold,                        
                                        const int *init_flag,
                                        const int *is_ghost);

/* ------------------------------- gradient threshold --------------------------------- */
#define FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_gradient_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH)

int FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH(const int *blockno,
                                          const double *qval, 
                                          const double* qmin, 
                                          const double *qmax,
                                          const double quad[], 
                                          const double *dx, 
                                          const double *dy, 
                                          const double *xc, 
                                          const double *yc, 
                                          const double *tag_threshold,
                                          const int *init_flag,
                                          const int *is_ghost);


/* ------------------------------- user threshold --------------------------------- */
#define USER_EXCEEDS_TH \
                  FCLAW_F77_FUNC(user_exceeds_th, \
                                 USER_EXCEEDS_TH)

int USER_EXCEEDS_TH(const int *blockno,
                    const double *qval, 
                    const double* qmin, 
                    const double *qmax,
                    const double quad[], 
                    const double *dx, 
                    const double *dy, 
                    const double *xc, 
                    const double *yc, 
                    const double *tag_threshold,
                    const int *init_flag,
                    const int *is_ghost);

    
#ifdef __cplusplus
}
#endif

#endif
