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

#include <fclaw2d_patch.hpp>
#include <fclaw2d_farraybox.hpp>


class ClawPatch
{
public :

    ClawPatch();
    ~ClawPatch();

    void define(fclaw2d_domain_t* domain,
                fclaw2d_patch_t* this_patch,
                int a_blockno,
                fclaw2d_build_mode_t build_mode);

    void copyFrom(ClawPatch *a_cp);

    // ----------------------------------------------------------------
    // Time stepping routines
    // ----------------------------------------------------------------

    void save_step();
    void save_current_step();
    void restore_step();

    void ghost_comm(double *qpack, int time_interp, int packmode);
    void partition_pack(double* qdata);
    void partition_unpack(double *qdata);
    void initialize_after_partition();
    void initialize_after_regrid();

    void setup_area_storage();
    void setup_metric_storage();

    void setup_for_time_interpolation(const double& alpha, const int& psize);
    // void reset_after_time_interpolation();
    double* q_time_sync(fclaw_bool time_interp);
    void finegrid_neighbors(int y);
    int has_finegrid_neighbors();
    FArrayBox newGrid();

    Box dataBox();  /* Box containing data pointer q */
    Box areaBox();  /* Box containing area */
    Box edgeBox();  /* Box containing edge based values */
    Box nodeBox();  /* Box containing nodes */

    /* ----------------------------------------------------------------
       Internal boundary conditions
       ---------------------------------------------------------------- */
    void exchange_face_ghost(const int& a_iface,
                             ClawPatch *a_neighbor_cp,
                             int time_interp,
                             fclaw2d_transform_data_t* transform_data);

    // Exchange corner ghost cells with interior neighbor (interior to the domain)
    void exchange_corner_ghost(const int& a_corner, ClawPatch *cp_neighbor,
                               int time_interp,
                               fclaw2d_transform_data_t* transform_data);

    void mb_exchange_block_corner_ghost(const int& a_icorner,
                                        ClawPatch *a_neighbor_cp,
                                        int time_interp);


    // Average finer grid values onto coarser grid ghost cells
    void average_face_ghost(const int& a_idir,
                            const int& a_iside,
                            const int& a_num_neighbors,
                            const int& a_refratio,
                            ClawPatch *neighbor_cp,
                            fclaw_bool a_time_interp,
                            const int& igrid,
                            fclaw2d_transform_data_t* transform_data);

    void average_corner_ghost(const int& a_corner, const int& a_refratio,
                              ClawPatch *cp_fine,
                              fclaw_bool a_time_interp,
                              fclaw2d_transform_data_t* transform_cptr);

    void mb_average_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                       ClawPatch *cp_fine,
                                       fclaw_bool a_time_interp);

    // Interpoalte coarser values onto fine grid ghost
    void interpolate_face_ghost(const int& a_idir,
                                const int& a_iside,
                                const int& a_num_neighbors,
                                const int& a_refratio,
                                ClawPatch *a_neighbor_cp,
                                fclaw_bool a_time_interp,
                                const int& igrid,
                                fclaw2d_transform_data_t* transform_data);

    void interpolate_corner_ghost(const int& a_corner, const int& a_refratio,
                                  ClawPatch *cp_fine,
                                  fclaw_bool a_time_interp,
                                  fclaw2d_transform_data_t* transform_cptr);

    void mb_interpolate_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                           ClawPatch *cp_fine,
                                           fclaw_bool a_time_interp);

    void set_boundary_to_value(const int& time_interp,
                               double& value);

    void set_corners_to_value(const int& time_interp,
                              double& value);

    /* ----------------------------------------------------------------
       Pillow grid ghost exchanges
       ---------------------------------------------------------------- */
    // ----------------------------------------------------------------------------------
    // Physical boundary conditions
    // ----------------------------------------------------------------------------------

    // Boundary condition routines that don't require neighbors
    void set_phys_face_ghost(const fclaw_bool a_intersects_bc[],
                             const int a_mthbc[],
                             const double& t,
                             const double& dt);

    // ----------------------------------------------------------------
    // Mapped grids
    // ----------------------------------------------------------------

    void setup_manifold(const int& a_level,
                        const amr_options_t *gparms,
                        fclaw2d_build_mode_t build_mode);

    void set_block_corner_count(const int icorner, const int block_corner_count);

    // ----------------------------------------------------------------
    // Miscellaneous
    // ----------------------------------------------------------------

    int size();

    // ----------------------------------------------------------------
    // Access functions
    // ----------------------------------------------------------------
    double mx();
    double my();
    double mbc();
    double meqn();

    double dx();
    double dy();

    double *xp();
    double *yp();
    double *zp();
    double *xd();
    double *yd();
    double *zd();

    double *area();

    double *xface_normals();
    double *yface_normals();
    double *xface_tangents();
    double *yface_tangents();
    double *surf_normals();
    double *curvature();
    double *edge_lengths();

    double xlower();
    double ylower();
    double xupper();
    double yupper();

    double* q();

    int* block_corner_count();

    void* clawpack_patch_data(int id);

    static fclaw_app_t* app;
    static fclaw2d_global_t *global;
    static int pack_layers;
    static int ghost_pack_size;
    static int ghost_patch_pack_area;

protected :

    /* Solution data */
    int m_meqn;                    /* also in amr_options_t */
    FArrayBox m_griddata;
    FArrayBox m_griddata_last;
    FArrayBox m_griddata_save;
    FArrayBox m_griddata_time_interpolated;

    /* Grid info */
    int m_mx;           /* also in amr_options_t */
    int m_my;           /* also in amr_options_t */
    int m_mbc;          /* also in amr_options_t */

    double m_dx;
    double m_dy;
    double m_xlower;
    double m_ylower;
    double m_xupper;
    double m_yupper;

    fclaw_bool m_manifold;    /* also in amr_options_t */
    int m_blockno;

    FArrayBox m_xp;
    FArrayBox m_yp;
    FArrayBox m_zp;

    FArrayBox m_xd;
    FArrayBox m_yd;
    FArrayBox m_zd;

    FArrayBox m_xface_normals;
    FArrayBox m_yface_normals;
    FArrayBox m_xface_tangents;
    FArrayBox m_yface_tangents;
    FArrayBox m_surf_normals;
    FArrayBox m_edge_lengths;

    FArrayBox m_area;
    FArrayBox m_curvature;  // ???

    int m_block_corner_count[4];
    int m_has_finegrid_neighbors;

    /* This is an opaque pointer */
    fclaw_package_data_t *m_package_data_ptr;

};

ClawPatch* fclaw2d_clawpatch_get_cp(fclaw2d_patch_t* this_patch);

#endif /* !FCLAW2D_CLAWPATCH_HPP */
