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

#ifndef __CLAWPATCH_H
#define __CLAWPATCH_H

#include "fclaw2d_farraybox.H"
#include "fclaw2d_clawpatch.H"
#include "fclaw2d_transform.h"
#include "fclaw_options.h"
#include <fclaw_package.h>

class ClawPatch
{
public :

    ClawPatch();
    ~ClawPatch();

    void define(const double&  a_xlower,
                const double&  a_ylower,
                const double&  a_xupper,
                const double&  a_yupper,
                const int& a_blockno,
                const int& a_level,
                const amr_options_t *a_gparms);

    void copyFrom(ClawPatch *a_cp);

    // ----------------------------------------------------------------
    // Time stepping routines
    // ----------------------------------------------------------------

    void save_step();
    void save_current_step();
    void restore_step();

    void pack_griddata(double* qdata);
    void unpack_griddata(double* qdata);
    void unpack_griddata_time_interpolated(double* qdata);

    void setup_for_time_interpolation(const double& alpha);
    void reset_after_time_interpolation();
    double* q_time_sync(fclaw_bool time_interp);
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
                             fclaw2d_transform_data_t* transform_data);

    // Average finer grid values onto coarser grid ghost cells
    void average_face_ghost(const int& a_idir,
                            const int& a_iside,
                            const int& a_num_neighbors,
                            const int& a_refratio,
                            ClawPatch *neighbor_cp,
                            fclaw_bool a_time_interp,
                            const int& igrid,
                            fclaw2d_transform_data_t* transform_data);

    // Interpoalte coarser values onto fine grid ghost
    void interpolate_face_ghost(const int& a_idir,
                                const int& a_iside,
                                const int& a_num_neighbors,
                                const int& a_refratio,
                                ClawPatch *a_neighbor_cp,
                                fclaw_bool a_time_interp,
                                const int& igrid,
                                fclaw2d_transform_data_t* transform_data);

    // Exchange corner ghost cells with interior neighbor (interior to the domain)
    void exchange_corner_ghost(const int& a_corner, ClawPatch *cp_neighbor,
                               fclaw2d_transform_data_t* transform_data);

    void average_corner_ghost(const int& a_corner, const int& a_refratio,
                              ClawPatch *cp_fine,
                              fclaw_bool a_time_interp,
                              fclaw2d_transform_data_t* transform_cptr);

    void interpolate_corner_ghost(const int& a_corner, const int& a_refratio,
                                  ClawPatch *cp_fine,
                                  fclaw_bool a_time_interp,
                                  fclaw2d_transform_data_t* transform_cptr);

    /* ----------------------------------------------------------------
       Pillow grid ghost exchanges
       ---------------------------------------------------------------- */
    void mb_exchange_block_corner_ghost(const int& a_icorner,
                                        ClawPatch *a_neighbor_cp);

    void mb_average_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                       ClawPatch *cp_fine,
                                       fclaw_bool a_time_interp);

    void mb_interpolate_block_corner_ghost(const int& a_corner, const int& a_refratio,
                                           ClawPatch *cp_fine,
                                           fclaw_bool a_time_interp);

    // ----------------------------------------------------------------------------------
    // Physical boundary conditions
    // ----------------------------------------------------------------------------------

    // Boundary condition routines that don't require neighbors
    void set_phys_face_ghost(const fclaw_bool a_intersects_bc[],
                             const int a_mthbc[],
                             const double& t,
                             const double& dt);

    // Corners that lie on a physical face - data has to be exchanged
    // in a particular way
    void set_phys_corner_ghost(const int& a_corner, const int a_mthbc[],
                               const double& t, const double& dt);

    // Exchange corner ghost with boundary neighbor
    void exchange_phys_face_corner_ghost(const int& a_corner,
                                         const int& a_iface, ClawPatch* cp);

    // ----------------------------------------------------------------
    // Mapped grids
    // ----------------------------------------------------------------

    void setup_manifold(const int& a_level, const amr_options_t *gparms);

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

    /* This is an opaque pointer */
    fclaw_package_data_t *m_package_data_ptr;

};

#endif
