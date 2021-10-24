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

#ifndef FCLAW2D_CLAWPATCH_TRANSFORM_H
#define FCLAW2D_CLAWPATCH_TRANSFORM_H

/**
 * @file 
 * Transforms clawpatch indexes
 */

#include <fclaw_base.h>   /* Defines FCLAW_F77_FUNC */


#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;  /* fclaw2d_patch.h includes this file */
/** Struture to store a patch's transform data */
struct fclaw2d_patch_transform_data;

/**
 * @brief Initialize the patch's tranform data
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in,out] transform the transform data
 */
void fclaw2d_clawpatch_transform_init_data(struct fclaw2d_global* glob, 
                                           struct fclaw2d_patch* this_patch,
                                           int blockno, int patchno,
                                           struct fclaw2d_patch_transform_data* transform);

/** @copydoc fclaw2d_patch_face_transformation() */
void fclaw2d_clawpatch_face_transformation (int faceno, int rfaceno, int ftransform[]);

/** @copydoc fclaw2d_patch_face_transformation_intra() */
void fclaw2d_clawpatch_face_transformation_intra (int ftransform[]);


/* ----------------------------- Called from Fortran ---------------------------------- */

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH_TRANSFORM_FACE FCLAW_F77_FUNC_(fclaw2d_clawpatch_transform_face, \
                                                         FCLAW2D_CLAWPATCH_TRANSFORM_FACE)
/**
 * @brief Tranform an index for a face-neighboring patch's coordinate system
 * 
 * @param[in] i1, j1 the input index 
 * @param[out] i2, j2 the transformed equivalent index in the neighboring patch
 * @param[in] ptdata the transform data
 */
void FCLAW2D_CLAWPATCH_TRANSFORM_FACE (const int *i1, const int *j1,
                                       int *i2, int *j2,
                                       struct fclaw2d_patch_transform_data** ptdata);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH_TRANSFORM_FACE_HALF FCLAW_F77_FUNC_(fclaw2d_clawpatch_transform_face_half, \
                                                              FCLAW2D_CLAWPATCH_TRANSFORM_FACE_HALF)
/**
 * @brief Tranform an index for a face-neighboring finer patch's coordinate system
 * 
 * @param[in] i1, j1 the input index 
 * @param[out] i2, j2 the four equivalent indexes in the finer patch
 * @param[in] ptdata the transform data
 */
void FCLAW2D_CLAWPATCH_TRANSFORM_FACE_HALF (const int *i1, const int *j1,
                                            int i2[], int j2[],
                                            struct fclaw2d_patch_transform_data** ptdata);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH_TRANSFORM_CORNER FCLAW_F77_FUNC_(fclaw2d_clawpatch_transform_corner, \
                                                           FCLAW2D_CLAWPATCH_TRANSFORM_CORNER)
/**
 * @brief Tranform an index for a face-neighboring patch's coordinate system
 * 
 * @param[in] i1, j1 the input index 
 * @param[out] i2, j2 the transformed equivalent index in the neighboring patch
 * @param[in] ptdata the transform data
 */
void FCLAW2D_CLAWPATCH_TRANSFORM_CORNER (const int *i1, const int *j1,
                                         int *i2, int *j2,
                                         struct fclaw2d_patch_transform_data** ptdata);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH_TRANSFORM_CORNER_HALF FCLAW_F77_FUNC_(fclaw2d_clawpatch_transform_corner_half, \
                                                                FCLAW2D_CLAWPATCH_TRANSFORM_CORNER_HALF)

/**
 * @brief Tranform an index for a corner-neighboring finer patch's coordinate system
 * 
 * @param[in] i1, j1 the input index 
 * @param[out] i2, j2 the four equivalent indexes in the finer patch
 * @param[in] ptdata the transform data
 */
void FCLAW2D_CLAWPATCH_TRANSFORM_CORNER_HALF (const int *i1, const int *j1,
                                              int i2[], int j2[],
                                              struct fclaw2d_patch_transform_data** ptdata);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_TRANSFORM_H */
