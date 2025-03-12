/*
Copyright (c) 2012-2025 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW_CLAWPATCH_GAUGES_H
#define FCLAW_CLAWPATCH_GAUGES_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw_global;
struct fclaw_gauge;
struct fclaw_patch;
struct fclaw_block;

/*  
    USER : This struct can be used to store any user defined information that 
    any gauge should have access to.  Right now, we store only q and aux values.

    Instances of this struct are stored as void* in the fclaw_gauge_t struct.
*/    

typedef struct fclaw_clawpatch_gauge_data
{
    int level;
    double tcurr;
    int meqn;      /* Needed when printing out buffers */
    int maux;
    double *qvar;  /* Store qvalues */
    double *avar;  /* Store aux variables */

} fclaw_clawpatch_gauge_data_t;


void fclaw_clawpatch_gauges_read_data(struct fclaw_global *glob, 
                                      struct fclaw_gauge **gauges, 
                                      int *num, int* dim);

void fclaw_clawpatch_gauges_create_files(struct fclaw_global *glob, 
                                         struct fclaw_gauge *gauges, 
                                         int num_gauges);

void fclaw_clawpatch_gauges_normalize_coordinates(struct fclaw_global *glob, 
                                                  struct fclaw_block *block,
                                                  int blockno, 
                                                  struct fclaw_gauge *g,
                                                  double *xc, double *yc, double *zc);

void fclaw_clawpatch_gauges_update(struct fclaw_global* glob, 
                                   struct fclaw_block* block,
                                   struct fclaw_patch* patch, 
                                   int blockno, int patchno,
                                   double tcurr, 
                                   struct fclaw_gauge *g);

void fclaw_clawpatch_gauges_print(struct fclaw_global *glob, 
                                  struct fclaw_gauge *gauge);

/* ----------------------------------- FORTRAN  ---------------------------------------- */

#define FCLAW2D_CLAWPATCH46_FORT_GAUGES_UPDATE \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_gauges_update, \
                          FCLAW2D_CLAWPATCH46_FORT_GAUGES_UPDATE)

void FCLAW2D_CLAWPATCH46_FORT_GAUGES_UPDATE(const int *num, 
                                           const int *mx, 
                                           const int *my, 
                                           const int* mbc, 
                                           const int *meqn, 
                                           const double* xlower,
                                           const double *ylower,
                                           const double *dx, 
                                           const double* dy,
                                           double q[],
                                           const int *maux, 
                                           double aux[],
                                           const double *xc, 
                                           double *yc,
                                           double qvar[], 
                                           double avar[]);


#define FCLAW2D_CLAWPATCH5_FORT_GAUGES_UPDATE \
           FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_gauges_update, \
                          FCLAW2D_CLAWPATCH5_FORT_GAUGES_UPDATE)
void FCLAW2D_CLAWPATCH5_FORT_GAUGES_UPDATE(const int *num, 
                                          const int *mx, 
                                          const int *my, 
                                          const int* mbc, 
                                          const int *meqn, 
                                          const double* xlower,
                                          const double *ylower,
                                          const double *dx, 
                                          const double* dy,
                                          double q[],
                                          const int *maux, 
                                          double aux[],
                                          const double *xc, 
                                          double *yc,
                                          double qvar[], 
                                          double avar[]);



#define FCLAW3D_CLAWPATCH46_FORT_GAUGES_UPDATE \
    FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_gauges_update, \
                    FCLAW3D_CLAWPATCH46_FORT_GAUGES_UPDATE)
void FCLAW3D_CLAWPATCH46_FORT_GAUGES_UPDATE(const int *num, 
                                           const int *mx, 
                                           const int *my, 
                                           const int *mz, 
                                           const int *mbc, 
                                           const int *meqn,
                                           const double *xlower,
                                           const double *ylower, 
                                           const double *zlower,
                                           const double *dx, 
                                           const double *dy, 
                                           const double *dz,
                                           double q[],
                                           const int *maux,
                                           double aux[],
                                           const double *xc,
                                           const double *yc, 
                                           const double *zc,
                                           double qvar[], 
                                           double avar[]);

#ifdef __cplusplus
}
#endif

#endif
