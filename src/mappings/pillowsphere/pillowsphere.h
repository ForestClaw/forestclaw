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

#ifndef PILLOWSPHERE_H
#define PILLOWSPHERE_H

#include <fclaw_base.h>    /* Needed for def. of FCLAW_F77_FUNC */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_glob;
struct fclaw2d_patch;
struct fclaw2d_patch_transform;


void pillow_copy_block_corner(struct fclaw2d_global* glob,
                              struct fclaw2d_patch* this_patch, 
                              struct fclaw2d_patch *corner_patch,
                              int icorner,
                              int time_interp,
                              struct fclaw2d_patch_transform_data 
                              *transform_data);

void pillow_average_block_corner(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch, 
                                 struct fclaw2d_patch *corner_patch,
                                 int icorner,
                                 int time_interp,
                                 struct fclaw2d_patch_transform_data 
                                 *transform_data);

void pillow_interpolate_block_corner(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* this_patch, 
                                     struct fclaw2d_patch *corner_patch,
                                     int icorner,
                                     int time_interp,
                                     struct fclaw2d_patch_transform_data 
                                     *transform_data);


/* --------------------------------- Fortran code ------------------------------------- */

#define FCLAW2D_FORT_MB_EXCHANGE_BLOCK_CORNER_GHOST \
               FCLAW_F77_FUNC(fclaw2d_fort_mb_exchange_block_corner_ghost, \
                              FCLAW2D_FORT_MB_EXCHANGE_BLOCK_CORNER_GHOST)

void FCLAW2D_FORT_MB_EXCHANGE_BLOCK_CORNER_GHOST(int* mx, int* my,
                                                 int* mbc, int* meqn,
                                                 double qthis[], 
                                                 double qneighbor[], 
                                                 int* icorner,
                                                 int* iblock);

// Averaging at block boundaries between coarse and fine grids.
#define FCLAW2D_FORT_MB_AVERAGE_BLOCK_CORNER_GHOST \
          FCLAW_F77_FUNC(fclaw2d_fort_mb_average_block_corner_ghost,\
                         FCLAW2D_FORT_MB_AVERAGE_BLOCK_CORNER_GHOST)

void  FCLAW2D_FORT_MB_AVERAGE_BLOCK_CORNER_GHOST(int* mx, int* my, int* mbc,
                                                 int* meqn, 
                                                 int* refratio, 
                                                 double qcoarse[],
                                                 double qfine[], 
                                                 double areacoarse[], 
                                                 double areafine[],
                                                 int* a_coarse_corner,
                                                 int* blockno);

// Averaging at block boundaries between coarse and fine grids.
#define FCLAW2D_FORT_MB_INTERPOLATE_BLOCK_CORNER_GHOST \
          FCLAW_F77_FUNC(fclaw2d_fort_mb_interpolate_block_corner_ghost, \
                         FCLAW2D_FORT_MB_INTERPOLATE_BLOCK_CORNER_GHOST)

void  FCLAW2D_FORT_MB_INTERPOLATE_BLOCK_CORNER_GHOST(int* mx, int* my, int* mbc,
                                                     int* meqn, int* refratio,
                                                     double qcoarse[],
                                                     double qfine[], 
                                                     int* icoarse_corner,
                                                     int* blockno);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif    