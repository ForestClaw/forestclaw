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

#ifndef FCLAW2D_DEFS_H
#define FCLAW2D_DEFS_H

/* Not nice for namespace cleanliness, but necessary to include mpi.h */
#include <sc.h>

/* Some convenient preprocessor directives *
 * TODO: Having static int in a .H file leads to multiple definitions.
 *       Rather use extern const int declares and define them in a .c file.
 */

#define FCLAW_SPACEDIM 2
static const int SpaceDim = FCLAW_SPACEDIM;

/* Number of faces to a patch. Changed from CUBEFACES to NUMFACES to
   avoid any confusion in the 2d case. */
#define FCLAW_NUMFACES (2 * FCLAW_SPACEDIM)
static const int NumFaces = FCLAW_NUMFACES;

#define FCLAW_P4EST_REFINE_FACTOR 2
static const int p4est_refineFactor = FCLAW_P4EST_REFINE_FACTOR;

#define FCLAW_NUM_CORNERS 4
static const int NumCorners = FCLAW_NUM_CORNERS;

#define FCLAW_NUM_SIBLINGS 4
static const int NumSiblings = FCLAW_NUM_SIBLINGS;

/* change to P4EST_BOOL ?
 * No, p4est does not know about bool and works with pure int */
#define fclaw_bool bool
#define fclaw_true true
#define fclaw_false false

#define fclaw_cptr unsigned long int

#define fclaw_nan NAN

#endif
