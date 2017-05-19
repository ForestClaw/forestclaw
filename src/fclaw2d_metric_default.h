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

#ifndef FCLAW2D_METRIC_DEFAULT_H
#define FCLAW2D_METRIC_DEFAULT_H

#include <fclaw2d_metric_default_fort.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

void fclaw2d_metric_compute_area(struct fclaw2d_global *glob,
                                 struct fclaw2d_patch* this_patch,
                                 int blockno,
                                 int patchno);

void fclaw2d_metric_area_set_ghost(struct fclaw2d_global *glob,
                                   struct fclaw2d_patch* this_patch,
                                   int blockno,
                                   int patchno);

void fclaw2d_metric_compute_area_exact(struct fclaw2d_global *glob,
                                       struct fclaw2d_patch *this_patch,
                                       int blockno,
                                       int patchno);

void fclaw2d_metric_area_set_ghost_exact(struct fclaw2d_global *glob,
                                         struct fclaw2d_patch* this_patch,
                                         int blockno,
                                         int patchno);

void fclaw2d_metric_setup_mesh(struct fclaw2d_global *glob,
                               struct fclaw2d_patch *this_patch,
                               int blockno,
                               int patchno);

/* Includes computation of normals and tangents to face, as well
   as surface normals and curvature */
void fclaw2d_metric_compute_normals(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch *this_patch,
                                    int blockno,
                                    int patchno);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
