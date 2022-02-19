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

#ifndef FCLAW2D_TO_3D_H
#define FCLAW2D_TO_3D_H

#ifdef FORESTCLAW2D_H
#error "The include files forestclaw2d.h and fclaw2d_to_3d.h cannot be combined"
#endif

#include <p4est_to_p8est.h>

/* redefine typedefs */
#define fclaw2d_patch_flags_t           fclaw3d_patch_flags_t
#define fclaw2d_patch_t                 fclaw3d_patch_t
#define fclaw2d_block_t                 fclaw3d_block_t
#define fclaw2d_domain_t                fclaw3d_domain_t

/* redefine enums */
#define FCLAW2D_PATCH_ON_BLOCK_FACE_0   FCLAW3D_PATCH_ON_BLOCK_FACE_0
#define FCLAW2D_PATCH_ON_BLOCK_FACE_1   FCLAW3D_PATCH_ON_BLOCK_FACE_1
#define FCLAW2D_PATCH_ON_BLOCK_FACE_2   FCLAW3D_PATCH_ON_BLOCK_FACE_2
#define FCLAW2D_PATCH_ON_BLOCK_FACE_3   FCLAW3D_PATCH_ON_BLOCK_FACE_3
#define FCLAW2D_PATCH_FIRST_SIBLING     FCLAW3D_PATCH_FIRST_SIBLING
#define FCLAW2D_PATCH_IS_GHOST          FCLAW3D_PATCH_IS_GHOST
#define FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY FCLAW3D_PATCH_ON_PARALLEL_BOUNDARY

/* redefine variables */
#define fclaw2d_patch_block_face_flags  fclaw3d_patch_block_face_flags
#define fclaw2d_smallest_h              fclaw3d_smallest_h

/* redefine functions */
#define fclaw2d_domain_global_maximum   fclaw3d_domain_global_maximum
#define fclaw2d_domain_global_sum       fclaw3d_domain_global_sum
#define fclaw2d_domain_barrier          fclaw3d_domain_barrier
#define fclaw2d_domain_dimension        fclaw3d_domain_dimension
#define fclaw2d_check_initial_level     fclaw3d_check_initial_level
#define fclaw2d_domain_new_unitsquare   fclaw3d_domain_new_unitcube

/* translations not found in p4est */
#ifndef p4est_wrap_new_unitsquare
#define p4est_wrap_new_unitsquare       p8est_wrap_new_unitcube
#endif

#endif /* !FCLAW2D_TO_3D_H */
