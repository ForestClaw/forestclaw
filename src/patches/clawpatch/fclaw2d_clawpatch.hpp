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

#ifndef FCLAW2D_CLAWPATCH_HPP
#define FCLAW2D_CLAWPATCH_HPP

#include <fclaw2d_farraybox.hpp>

struct fclaw2d_patch;
struct fclaw2d_transform_data;

class fclaw2d_clawpatch_t
{
public :
    Box dataBox();  /* Box containing data pointer q */
    Box areaBox();  /* Box containing area */
    Box edgeBox();  /* Box containing edge based values */
    Box nodeBox();  /* Box containing nodes */

    /* ----------------------------------------------------------------
       Pillow grid ghost exchanges
       ---------------------------------------------------------------- */

    void mb_exchange_block_corner_ghost(const int& a_icorner,
                                        fclaw2d_clawpatch_t *a_neighbor_cp,
                                        int time_interp);

    void mb_average_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                       fclaw2d_clawpatch_t *cp_fine,
                                       int a_time_interp);

    void mb_interpolate_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                           fclaw2d_clawpatch_t *cp_fine,
                                           int a_time_interp);

    // ----------------------------------------------------------------
    // Miscellaneous
    // ----------------------------------------------------------------

    void* clawpack_patch_data(int id);

#if 0
    static fclaw_app_t* app;
    static fclaw2d_global_t *global;
#endif    

    /* Solution data */
    int meqn;                    /* also in amr_options_t */
    FArrayBox griddata;
    FArrayBox griddata_last;
    FArrayBox griddata_save;
    FArrayBox griddata_time_interpolated;
    FArrayBox griderror;

    /* Grid info */
    int mx;           /* also in amr_options_t */
    int my;           /* also in amr_options_t */
    int mbc;          /* also in amr_options_t */
    int maux;

    double dx;
    double dy;
    double xlower;
    double ylower;
    double xupper;
    double yupper;

    int manifold;    /* also in amr_options_t */
    int blockno;

    FArrayBox aux;

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

fclaw2d_clawpatch_t* fclaw2d_clawpatch_get_cp(struct fclaw2d_patch* this_patch);

#endif /* !FCLAW2D_CLAWPATCH_HPP */
