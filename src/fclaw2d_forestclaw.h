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

#ifndef FCLAW2D_FORESTCLAW_H
#define FCLAW2D_FORESTCLAW_H

#include <fclaw_forestclaw.h>

/* Basic objects */
#include "fclaw2d_global.h"


#include "forestclaw2d.h"
#include "fclaw2d_domain.h"
#include "fclaw2d_block.h"
#include "fclaw2d_patch.h"


/* Basic header files that are probably required by all apps */
#include <fclaw2d_options.h>
#include "fclaw2d_vtable.h"

/* Mapping interface - needed even if 'nomap' is used.  */
#include "p4est_connectivity.h"
#include "fclaw2d_map.h"
#include "fclaw2d_map_query.h"

#include <fclaw2d_diagnostics.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* Some convenient preprocessor directives */

#define FCLAW2D_SPACEDIM 2
extern const int fclaw2d_SpaceDim;

#define FCLAW2D_NUMFACES (2 * FCLAW2D_SPACEDIM)
extern const int fclaw2d_NumFaces;

#define FCLAW2D_REFINE_FACTOR 2
extern const int fclaw2d_RefineFactor;

#define FCLAW2D_NUM_CORNERS 4
extern const int fclaw2d_NumCorners;

#define FCLAW2D_NUM_SIBLINGS 4
extern const int fclaw2d_NumSiblings;

void fclaw2d_initialize(fclaw2d_global_t *glob);
void fclaw2d_run(fclaw2d_global_t *glob);
void fclaw2d_finalize(fclaw2d_global_t *glob);

void fclaw2d_after_regrid(fclaw2d_global_t *glob);
void fclaw2d_problem_setup(fclaw2d_global_t *glob);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif


/* Note: either we make this a C .h file, or we remove the extern "C". */

#endif
