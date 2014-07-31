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

#ifndef CLAWPACK_FORT_H
#define CLAWPACK_FORT_H

/* this header file must come first */
#include "fclaw2d_defs.H"

#include "fclaw2d_convenience.h"
#include "forestclaw2d.h"
#include "amr_options.h"

#include <iostream>
#include <cstdlib>
#include <vector>

#ifdef __cplusplus
extern "C"
{          /* Beginning of extern "C" */
#if 0
}
#endif
#endif

void set_common_levels_(const int& maxlevel,const int& a_level,
                        const int& refratio);

/* ----------------------------------------------------------------------------------
   Initialization routines
   ---------------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------ */
/* This is here for now;  may be taken out once I figure out Carsten's
   'mapping context' idea */
void set_maptype_cart_();
void set_maptype_sphere_();
void set_maptype_hemisphere_();
void set_maptype_disk_();
void set_map_defaults_();

int isflat_();
int iscart_();
int isdisk_();
int ishemisphere_();
int issphere_();
/* ------------------------------------------------------------------------ */

void set_block_(const int * a_blockno);

void average_face_ghost_mapped_(const int& mx, const int& my, const int& mbc,
                                const int& meqn, double qfine[],double qcoarse[],
                                const double areacoarse[], const double areafine[],
                                const int& idir, const int& iside,
                                const int& num_neighbors, const int& refratio,
                                const int& igrid);

void mb_average_face_ghost_(const int& mx, const int& my, const int& mbc,
                            const int& meqn, double qfine[],double qcoarse[],
                            const double areacoarse[], const double areafine[],
                            const int& idir, const int& iside, const int& num_neighbors,
                            const int& refratio, const int& igrid,
                            fclaw_cptr transform_cptr);

void fixcapaq2_(const int& mx, const int& my, const int& mbc, const int& meqn,
                double qcoarse[], double qfine[], double areacoarse[], double areafine[],
                const int& p4est_refineFactor,
                const int& refratio, const int& igrid);

/* ----------------------------------------------------------------------------------
   Internal boundary conditions
   ---------------------------------------------------------------------------------- */
void exchange_face_ghost_(const int& mx, const int& my, const int& mbc, const int& meqn,
                          double qthis[],double qneighbor[], const int& a_idir,
                          fclaw_cptr transform_cptr);

void average_face_ghost_(const int& mx, const int& my, const int& mbc,
                         const int& meqn,
                         double qcoarse[],double qfine[],
                         double areacoarse[], double areafine[],
                         const int& idir, const int& iside,
                         const int& num_neighbors,
                         const int& refratio, const int& igrid,
                         const int& manifold, fclaw_cptr transform_cptr);

void interpolate_face_ghost_(const int& mx, const int& my, const int& mbc,
                             const int& meqn,
                             double qthis[],double qcoarse[],
                             const int& idir, const int& iside,
                             const int& num_neighbors,
                             const int& refratio, const int& igrid,
                             fclaw_cptr transform_cptr);


void exchange_corner_ghost_(const int& mx, const int& my, const int& mbc,
                            const int& meqn, double this_q[],double neighbor_q[],
                            const int& a_corner,fclaw_cptr transform_cptr);

void average_corner_ghost_(const int& mx, const int& my, const int& mbc,
                           const int& meqn, const int& a_refratio,
                           double qcoarse[], double qfine[],
                           double areacoarse[], double areafine[],
                           const int& manifold,
                           const int& a_corner, fclaw_cptr transform_cptr);

void interpolate_corner_ghost_(const int& mx, const int& my, const int& mbc,
                               const int& meqn, const int& a_refratio, double this_q[],
                               double neighbor_q[], const int& a_corner,
                               fclaw_cptr transform_cptr);

void mb_exchange_face_ghost_(const int& mx, const int& my, const int& mbc, const int& meqn,
                             double qthis[], double qneighbor[], const int& iface,
                             const int& iblock);

void mb_exchange_corner_ghost_(const int& mx, const int& my, const int& mbc,
                               const int& meqn,
                               double qthis[], double qneighbor[], const int& icorner,
                               int bdry[], const int& iblock);

void mb_exchange_block_corner_ghost_(const int& mx, const int& my,
                                     const int& mbc, const int& meqn,
                                     double qthis[], double qneighbor[], const int& icorner,
                                     const int& iblock);

void mb_interpolate_face_ghost_(const int& mx, const int& my, const int& mbc,
                                const int& meqn,
                                double qthis[],double qcoarse[], const int& idir,
                                const int& iside,
                                const int& num_neighbors,const int& refratio,
                                const int& igrid);

void mb_average_corner_ghost_(const int& mx, const int& my, const int& mbc,
                              const int& meqn,
                              const int& refratio, double qcoarse[],  double qfine[],
                              double areacoarse[], double areafine[],
                              const int& icorner, int intersects_block[]);

// Averaging at block boundaries between coarse and fine grids.
void  mb_average_block_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                     const int& meqn,const int& refratio, double qcoarse[],
                                     double qfine[],double areacoarse[], double areafine[],
                                     const int& a_coarse_corner,
                                     const int& blockno);


void mb_interpolate_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                  const int& meqn, const int& refratio,
                                  double qcoarse[],  double qfine[], const int& icorner,
                                  int intersects_block[]);

// Averaging at block boundaries between coarse and fine grids.
void  mb_interpolate_block_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                         const int& meqn,const int& refratio,
                                         double qcoarse[],
                                         double qfine[],const int& a_coarse_corner,
                                         const int& blockno);


/* ----------------------------------------------------------------------------------
   Physical boundary conditions
   ---------------------------------------------------------------------------------- */
void set_phys_corner_ghost_(const int& mx, const int& my, const int& mbc,
                            const int& meqn, double q[],const int& icorner,
                            const double &t, const double& dt, const int mthbc[]);

void exchange_phys_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                 const int& meqn, double qthis[],double qneighbor[],
                                 const int& icorner, const int& iside);

/* ----------------------------------------------------------------------------------
   Tagging for refinement/coarsening
   ---------------------------------------------------------------------------------- */
void interpolate_to_fine_patch_(const int& mx,const int& my,const int& mbc,
                                const int& meqn, double qcoarse[], double qfine[],
                                const int& p4est_refineFactor,
                                const int& refratio, const int& igrid);

void average_to_coarse_patch_(const int& mx,const int& my,const int& mbc,
                              const int& meqn,
                              double qcoarse[],double qfine[],
                              double areacoarse[],double areafine[],
                              const int& p4est_refineFactor,
                              const int& refratio, const int& igrid,
                              const int& manifold);

#if 0
void average_to_coarse_mapped_(const int& mx,const int& my, const int& mbc,
                               const int& meqn, double qcoarse[], double qfine[],
                               const double areacoarse[], const double areafine[],
                               const int& p4est_refineFactor,
                               const int& refratio, const int& igrid);
#endif

/* ----------------------------------------------------------------------------------
   Mapped grids
   area                : 1 float
   xp,yp,zp, xd,yd zp  : 6 floats
   xnormals, ynormals  : 6 floats
   ytangents, ytangent : 6 floats
   edgelengths         : 2 floats
   surfnormals         : 3 floats
   curvature           : 1 float
   -----------------------------
   total               : 25 floats (in 2d!!)

   ---------------------------------------------------------------------------------- */
void setup_mesh_(const int& mx, const int& my, const int& mbc,
                 const double& xlower, const double& ylower,
                 const double& dx, const double& dy,
                 double xp[], double yp[], double zp[], double xd[], double yd[], double zd[]);

void compute_area_(const int& mx, const int& my, const int& mbc,
                   const double& m_dx, const double& m_dy,const double& m_xlower,
                   const double& m_ylower, double area[],
                   const int& level, const int& maxlevel, const int& refratio);

void compute_normals_(const int& mx, const int& my, const int& mbc,
                      double xp[], double yp[], double zp[],
                      double xd[], double yd[], double zd[],
                      double xnormals[],double ynormals[]);

void compute_tangents_(const int& mx, const int& my, const int& mbc,
                      double xd[], double yd[], double zd[],
                      double xtangents[],double ytangents[],double edge_lengths[]);

void compute_surf_normals_(const int& mx, const int& my, const int& mbc,
                           double xnormals[],double ynormals[],double edge_lengths[],
                           double curvature[], double surfnormals[], double area[]);


/* ----------------------------------------------------------------------------------
   Output and diagnostics
   ---------------------------------------------------------------------------------- */

void write_tfile_(const int &iframe,const double &time,const int &meqn,
                  const int &ngrids, const int& maux);

void new_qfile_(const int &iframe);

void write_qfile_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
                  const int& mx, const int& my, const double& xlower, const double& ylower,
                  const double& dx, const double& dy, double q[], const int& iframe,
                  const int& patch_idx,const int& level,const int& blockno);

void compute_sum_(const int& mx, const int& my, const int& mbc, const int& meqn,
                  const double& dx,
                  const double& dy, double q[], double& sum);

void interpdata_(const int& m,double qlast[],double qcurr[], double qinterp[],
                 const double& alpha);


int get_debug_flag_();
void set_debug_flag_(const int& iflag);
void set_debug_on_();
void set_debug_off_();
fclaw_bool debug_is_on_();
fclaw_bool debug_is_off_();

void set_debug_info_(const int& block_idx, const int& patch_idx, const int& level);
void reset_debug_info_();
int get_block_idx_();
int get_patch_idx_();
int get_level_();
void debug_here_(const int& n, const char str[]);

void dump_patch_(const int& mx, const int& my, const int& m_mbc, const int& meqn,
                 const int& mq, double q[]);

#ifdef __cplusplus
#if 0
{
#endif
}           /* end of extern "C" */
#endif


class Box
{
public:
    Box();
    Box(const int ll[], const int ur[]);
    Box(const Box& a_box);
    int smallEnd(int idir) const;
    int bigEnd(int idir) const;

private:
    int m_ll[2];
    int m_ur[2];
};

class FArrayBox
{
public:

    FArrayBox();
    FArrayBox(const FArrayBox& A);
    ~FArrayBox();
    void define(const Box& a_box, int a_fields);
    double* dataPtr();
    Box box();
    int fields();
    int size();
    void operator=(const FArrayBox& fbox);
    void copyToMemory(double *data);
    void copyFromMemory(double *data);
private:
    double *m_data;
    int m_size;
    Box m_box;
    int m_fields;
    void set_dataPtr(int size);

};

#endif
