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

#include <fclaw2d_global.h>
#include <fclaw2d_forestclaw.h>
#include <fclaw_options.h>
#include <fclaw_package.h>
#include <fclaw2d_transform.h>
#include <fclaw2d_clawpatch.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_farraybox.hpp>

class fclaw2d_clawpatch_t;


class fclaw2d_clawpatch_t
{
public :

    fclaw2d_clawpatch_t();
    ~fclaw2d_clawpatch_t();

    void define(fclaw2d_domain_t* domain,
                fclaw2d_patch_t* this_patch,
                int a_blockno,
                fclaw2d_build_mode_t build_mode);

#if 0
    void copyFrom(fclaw2d_clawpatch_t *a_cp);
#endif

    void save_step();
    void save_current_step();
    void restore_step();

    void ghost_comm(double *qpack, int time_interp, int packmode);
    void partition_pack(double* qdata);
    void partition_unpack(double *qdata);

    void setup_area_storage();
    void setup_metric_storage();

    // void reset_after_time_interpolation();
    double* q_time_sync(fclaw_bool time_interp);

    FArrayBox newGrid();

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
                                       fclaw_bool a_time_interp);

    void mb_interpolate_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                           fclaw2d_clawpatch_t *cp_fine,
                                           fclaw_bool a_time_interp);

    // ----------------------------------------------------------------
    // Mapped grids
    // ----------------------------------------------------------------

    void setup_manifold(const int& a_level,
                        const amr_options_t *gparms,
                        fclaw2d_build_mode_t build_mode);

    // ----------------------------------------------------------------
    // Miscellaneous
    // ----------------------------------------------------------------

    int size();

    // ----------------------------------------------------------------
    // Access functions
    // ----------------------------------------------------------------
    // double mx();
    // double my();
    // double mbc();
    // double meqn();

    // double dx();
    // double dy();

    // double *xp();
    // double *yp();
    // double *zp();
    // double *xd();
    // double *yd();
    // double *zd();

    // double *area();

    // double *xface_normals();
    // double *yface_normals();
    // double *xface_tangents();
    // double *yface_tangents();
    // double *surf_normals();
    // double *curvature();
    // double *edge_lengths();

    // double xlower();
    // double ylower();
    // double xupper();
    // double yupper();

    // double* q();
    // double* q_last();
    // double* q_timeinterp();
    // double* error();

    void* clawpack_patch_data(int id);

    static fclaw_app_t* app;
    static fclaw2d_global_t *global;

    // fclaw2d_clawpatch_t *clawp;

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

    double dx;
    double dy;
    double xlower;
    double ylower;
    double xupper;
    double yupper;

    fclaw_bool manifold;    /* also in amr_options_t */
    int blockno;

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

    /* This is an opaque pointer */
    fclaw_package_data_t *package_data_ptr;

};
#if 0
struct fclaw2d_clawpatch_t
{
public:
  int meqn;

  int mx;
  int my;
  int mbc;

  double dx;
  double dy;
  double xlower;
  double ylower;
  double xupper;
  double yupper;

  int manifold;
  int blockno;

  FArrayBox griddata;

  FArrayBox xp;
  FArrayBox yp;
  FArrayBox zp;

  FArrayBox xd;
  FArrayBox yd;
  FArrayBox zd;
};
#endif
fclaw2d_clawpatch_t* fclaw2d_clawpatch_get_cp(fclaw2d_patch_t* this_patch);

#endif /* !FCLAW2D_CLAWPATCH_HPP */
