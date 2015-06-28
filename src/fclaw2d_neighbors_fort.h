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
#include <fclaw2d_defs.h>

#include <forestclaw2d.h>
#include <fclaw2d_convenience.h>

#include <fclaw_options.h>
#include <fclaw2d_transform.h>

#ifdef __cplusplus
extern "C"
{          /* Beginning of extern "C" */
#if 0
}
#endif
#endif

/* ------------------------------------------------------------------------ */


void mb_average_face_ghost_(const int& mx, const int& my, const int& mbc,
                            const int& meqn, double qfine[],double qcoarse[],
                            const double areacoarse[], const double areafine[],
                            const int& idir, const int& iside, const int& num_neighbors,
                            const int& refratio, const int& igrid,
                            fclaw2d_transform_data_t** transform_cptr);

/* ----------------------------------------------------------------------------------
   Internal boundary conditions
   ---------------------------------------------------------------------------------- */
void exchange_face_ghost_(const int& mx, const int& my, const int& mbc, const int& meqn,
                          double qthis[],double qneighbor[], const int& a_idir,
                          fclaw2d_transform_data_t** transform_cptr);

void average_face_ghost_(const int& mx, const int& my, const int& mbc,
                         const int& meqn,
                         double qcoarse[],double qfine[],
                         double areacoarse[], double areafine[],
                         const int& idir, const int& iside,
                         const int& num_neighbors,
                         const int& refratio, const int& igrid,
                         const int& manifold, fclaw2d_transform_data_t** transform_cptr);

void interpolate_face_ghost_(const int& mx, const int& my, const int& mbc,
                             const int& meqn,
                             double qthis[],double qcoarse[],
                             const int& idir, const int& iside,
                             const int& num_neighbors,
                             const int& refratio, const int& igrid,
                             fclaw2d_transform_data_t** transform_cptr);


void exchange_corner_ghost_(const int& mx, const int& my, const int& mbc,
                            const int& meqn, double this_q[],double neighbor_q[],
                            const int& a_corner,fclaw2d_transform_data_t** transform_cptr);

void average_corner_ghost_(const int& mx, const int& my, const int& mbc,
                           const int& meqn, const int& a_refratio,
                           double qcoarse[], double qfine[],
                           double areacoarse[], double areafine[],
                           const int& manifold,
                           const int& a_corner, fclaw2d_transform_data_t** transform_cptr);

void interpolate_corner_ghost_(const int& mx, const int& my, const int& mbc,
                               const int& meqn, const int& a_refratio, double this_q[],
                               double neighbor_q[], const int& a_corner,
                               fclaw2d_transform_data_t** transform_cptr);

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
   Mapped grids
   area                : 1 float
   xp,yp,zp, xd,yd zp  : 6 floats
   xnormals, ynormals  : 6 floats
   ytangents, ytangent : 6 floats
   edgelengths         : 2 floats
   surfnormals         : 3 floats
   curvature           : 1 float
   -----------------------------
   total               : 25 floats per field (in 2d!!)

   ---------------------------------------------------------------------------------- */
void setup_mesh_(const int& mx, const int& my, const int& mbc,
                 const double& xlower, const double& ylower,
                 const double& dx, const double& dy,int& blockno,
                 double xp[], double yp[], double zp[],
                 double xd[], double yd[], double zd[]);

void compute_area_(const int& mx, const int& my, const int& mbc,
                   const double& m_dx, const double& m_dy,const double& m_xlower,
                   const double& m_ylower, const int& blockno, double area[],
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

#ifdef __cplusplus
#if 0
{
#endif
}           /* end of extern "C" */
#endif

#endif
