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

#ifndef FCLAW2D_TRANSFORM_H
#define FCLAW2D_TRANSFORM_H


#include "fclaw2d_typedefs.h"
#include "forestclaw2d.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


/* Data structure for transformation data at non-aligned boundaries */
typedef struct fclaw2d_transform_data
{
    fclaw2d_patch_t* this_patch;
    fclaw2d_patch_t* neighbor_patch;
    int transform[9];
    int mx;
    int my;
    /* This seems to be okay - we should never need both a corner and a face */
    int iface;
    int icorner;
    int based;   // 0 for cell-centered; 1 for nodes
    int fine_grid_pos;   // Fine grid position (0 or 1)
} fclaw2d_transform_data_t;


/* Underscore is needed since these routines will be called from Fortran */
void transform_face_samesize_(const int &i1, const int &j1,
                              int *i2, int *j2,
                              fclaw2d_transform_data_t* tdata);

void transform_face_halfsize_(const int &i1, const int &j1,
                              int i2[], int j2[],
                              fclaw2d_transform_data_t* tdata);


/* Underscore is needed since these routines will be called from Fortran */
void transform_corner_samesize_(const int &i1, const int &j1,
                                int *i2, int *j2,
                                fclaw2d_transform_data_t* tdata);

void transform_corner_halfsize_(const int &i1, const int &j1,
                                int i2[], int j2[],
                                fclaw2d_transform_data_t* tdata);

void transform_corner_halfsize2_(const int &i1, const int &j1,
                                int i2[], int j2[],
                                fclaw2d_transform_data_t* tdata);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
