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

#ifndef FCLAW2D_METRIC_HPP
#define FCLAW2D_METRIC_HPP

#include <fclaw2d_farraybox.hpp>     /* Needed for FArray box used to store metric data */

typedef struct fclaw2d_metric_vtable fclaw2d_metric_vtable_t;

struct fclaw2d_global;
struct fclaw2d_patch;

class fclaw2d_metric_patch_t
{
public :

    int mx;           
    int my;           
    int mbc;
    int blockno;

    double dx;
    double dy;
    double xlower;
    double ylower;
    double xupper;
    double yupper;
    
    FArrayBox xp;
    FArrayBox yp;
    FArrayBox zp;

    FArrayBox xd;
    FArrayBox yd;
    FArrayBox zd;

    FArrayBox xface_normals;
    FArrayBox yface_normals;
    FArrayBox xface_tangents;
    FArrayBox yface_tangents;
    FArrayBox surf_normals;
    FArrayBox edge_lengths;

    FArrayBox area;
    FArrayBox curvature;  // ???
};

fclaw2d_metric_patch_t* fclaw2d_metric_patch_new();

void fclaw2d_metric_patch_delete(fclaw2d_metric_patch_t **patchmp);

#endif /* !FCLAW2D_METRIC_HPP */
