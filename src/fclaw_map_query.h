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


#ifndef FCLAW_MAP_QUERY_H
#define FCLAW_MAP_QUERY_H

#include <fclaw_map.h>
#include <fclaw_map_query_defs.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

#define FCLAW_MAP_IS_USED FCLAW_F77_FUNC_(fclaw_map_is_used,FCLAW_MAP_IS_USED)
int FCLAW_MAP_IS_USED(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_CART FCLAW_F77_FUNC_(fclaw_map_is_cart,FCLAW_MAP_IS_CART)
int FCLAW_MAP_IS_CART(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_AFFINE FCLAW_F77_FUNC_(fclaw_map_is_affine,FCLAW_MAP_IS_AFFINE)
int FCLAW_MAP_IS_AFFINE(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_DISK FCLAW_F77_FUNC_(fclaw_map_is_disk,FCLAW_MAP_IS_DISK)
int FCLAW_MAP_IS_DISK(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_PILLOWDISK FCLAW_F77_FUNC_(fclaw_map_is_pillowdisk,FCLAW_MAP_IS_PILLOWDISK)
int FCLAW_MAP_IS_PILLOWDISK(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_SQUAREDDISK FCLAW_F77_FUNC_(fclaw_map_is_squareddisk,FCLAW_MAP_IS_SQUAREDDISK)
int FCLAW_MAP_IS_SQUAREDDISK(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_PILLOWSPHERE FCLAW_F77_FUNC_(fclaw_map_is_pillowsphere,FCLAW_MAP_IS_PILLOWSPHERE)
int FCLAW_MAP_IS_PILLOWSPHERE(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_CUBEDSPHERE FCLAW_F77_FUNC_(fclaw_map_is_cubedsphere,FCLAW_MAP_IS_CUBEDSPHERE)
int FCLAW_MAP_IS_CUBEDSPHERE(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_FLAT FCLAW_F77_FUNC_(fclaw_map_is_flat,FCLAW_MAP_IS_FLAT)
int FCLAW_MAP_IS_FLAT(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_SPHERE FCLAW_F77_FUNC_(fclaw_map_is_sphere,FCLAW_MAP_IS_SPHERE)
int FCLAW_MAP_IS_SPHERE(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_HEMISPHERE FCLAW_F77_FUNC_(fclaw_map_is_hemisphere,FCLAW_MAP_IS_HEMISPHERE)
int FCLAW_MAP_IS_HEMISPHERE(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_TORUS FCLAW_F77_FUNC_(fclaw_map_is_torus,FCLAW_MAP_IS_TORUS)
int FCLAW_MAP_IS_TORUS(fclaw_map_context_t** cont);

#define FCLAW_MAP_IS_BRICK FCLAW_F77_FUNC_(fclaw_map_is_brick,FCLAW_MAP_IS_BRICK)
int FCLAW_MAP_IS_BRICK(fclaw_map_context_t** cont);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif


#endif
