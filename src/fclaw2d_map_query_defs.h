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

#ifndef FCLAW2D_MAP_QUERY_DEFS_H
#define FCLAW2D_MAP_QUERY_DEFS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* These integers are meant to be passed in query_identifier of map_query_t.
   Note : These categories are not completely well defined, and are used
   primarily for the user's benefit and so subject to the user's interpretation
*/
#define FCLAW2D_MAP_QUERY_IS_USED           0     /* is the map used at all? */
#define FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT    1     /* x_i -> a_i x_i + b_i */
#define FCLAW2D_MAP_QUERY_IS_AFFINE         2     /* x -> A x + b */
#define FCLAW2D_MAP_QUERY_IS_NONLINEAR      3     /* x -> F(x) */
#define FCLAW2D_MAP_QUERY_IS_CART           4     /* x -> x   (Cartesian) */
#define FCLAW2D_MAP_QUERY_IS_GRAPH          5     /* (x,y) -> (x,y,f(x,y)) */
#define FCLAW2D_MAP_QUERY_IS_PLANAR         6     /* (x,y) -> (?,?,0) */
#define FCLAW2D_MAP_QUERY_IS_ALIGNED        7     /* (x,y) -> (f(x),g(y),0) */
#define FCLAW2D_MAP_QUERY_IS_FLAT           8     /* Zero curvature */
#define FCLAW2D_MAP_QUERY_IS_SPHERE         9     /* A topological sphere */
#define FCLAW2D_MAP_QUERY_IS_DISK          10     /* A topological disk */
#define FCLAW2D_MAP_QUERY_IS_PILLOWDISK    11     /* Pillow disk */
#define FCLAW2D_MAP_QUERY_IS_SQUAREDDISK   12     /* Squared disk */
#define FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE  13     /* Pillow sphere */
#define FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE   14     /* Cubed sphere */
#define FCLAW2D_MAP_QUERY_IS_FIVEPATCH     15     /* Five patch unit square */
#define FCLAW2D_MAP_QUERY_IS_BILINEAR      16     /* Five patch unit square */
#define FCLAW2D_MAP_QUERY_IS_HEMISPHERE    17     /* Hemisphere grid */
#define FCLAW2D_MAP_QUERY_IS_TORUS         18     /* Hemisphere grid */
#define FCLAW2D_MAP_QUERY_IS_BRICK         19     /* Is a Cartesian brick */
#define FCLAW2D_MAP_QUERY_LAST             20     /* Number of "official" queries. */

#if 0
/* Generic query function (kept up to date with list above) */
static int
fclaw2d_map_query_generic (fclaw2d_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW2D_MAP_QUERY_IS_USED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_AFFINE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_NONLINEAR:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PLANAR:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_FLAT:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_DISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PILLOWDISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SQUAREDDISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_FIVEPATCH:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_BILINEAR:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_HEMISPHERE:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_generic (fclaw2d_map_query_defs.h) : " \
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}
#endif

#if 0
/* Sample list used for Fortran map.  Note : It is not recommended that you
   use this form, as it is hard to keep lists used in applications consistent
   with list above.*/
  int query_results[FCLAW2D_MAP_QUERY_LAST];
  query_results[FCLAW2D_MAP_QUERY_IS_USED] = 1;
  query_results[FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_AFFINE] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_NONLINEAR] = 1;
  query_results[FCLAW2D_MAP_QUERY_IS_CART] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_GRAPH] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_PLANAR] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_ALIGNED] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_FLAT] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_SPHERE] = 1;
  query_results[FCLAW2D_MAP_QUERY_IS_DISK] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_PILLOWDISK] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_SQUAREDDISK] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_FIVEPATCH] = 0;
  query_results[FCLAW2D_MAP_QUERY_IS_HEMISPHERE] = 0;
#endif

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
