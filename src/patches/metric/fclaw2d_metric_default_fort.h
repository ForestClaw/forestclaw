/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW2D_METRIC_DEFAULT_FORT_H
#define FCLAW2D_METRIC_DEFAULT_FORT_H

#include "forestclaw2d.h"   /* Needed for FCLAW_F77_FUNC macro */

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


/**
 * @file
 * C declarations for Fortran subroutines
 */


/* Typedef'd functions that user can set to whatever they want */
/**
 * @brief Compute the cell center and node coordinates for a patch
 *
 * @param[in] mx, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] xlower, ylower the lower left coordinate of the patch
 * @param[in] dx, dy the spacings in the x and y direcitons
 * @param[in] blockno the block number
 * @param[out] xp, yp, zp the coordinates of cell centers
 * @param[out] xd, yd, zd the coordinates nodes
 */
typedef void (*fclaw2d_metric_fort_compute_mesh_t)(const int* mx, const int* my, 
                                                   const int* mbc,
                                                   const double* xlower, 
                                                   const double* ylower,
                                                   const double* dx, const double* dy,
                                                   int* blockno,
                                                   double xp[], double yp[], double zp[],
                                                   double xd[], double yd[], double zd[]);
/**
 * @brief Compute the area for each cell
 *
 * @param[in] mx, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] dx, dy the spacings in the x and y direcitons
 * @param[in] xlower, ylower the lower left coordinate of the patch
 * @param[in] blockno the block number
 * @param[out] area the area of each cell
 * @param[in] quadsize the length of the quad
 * @param[in] quadstore stores a group of cell values
 * @param[in] ghost_only true if only ghost cell shouyld be computed
 */
typedef void (*fclaw2d_metric_fort_compute_area_t)(const int* mx, 
                                                   const int* my, const int* mbc,
                                                   const double* dx, const double* dy,
                                                   const double* xlower, 
                                                   const double* ylower,
                                                   const int* blockno, double area[],
                                                   const int* level, const int* maxlevel,
                                                   const int* refratio, 
                                                   const int* ghost_only);

/**
 * @brief Compute the normals at the cell faces
 *
 * @param[in] mx, my the number of cells in the x and y direcitons
 * @param[in] mbc the number of ghost cells
 * @param[in] xp, yp, zp the cell center coordinates
 * @param[in] xd, yd, zd the node coordinates
 * @param[in] xnormals, ynormals the normals of the cell faces
 */
typedef void (*fclaw2d_metric_fort_compute_normals_t)(const int* mx, const int* my, 
                                                      const int* mbc,
                                                      double xp[], double yp[], 
                                                      double zp[],
                                                      double xd[], double yd[], 
                                                      double zd[],
                                                      double xnormals[], double ynormals[]);

/**
 * @brief Computes the tangents and lengths at the cell faces
 *
 * @param[in] mx, my the number of cells in the x and y direcitons
 * @param[in] mbc the number of ghost cells
 * @param[in] xd, yd, zd the node coordinates
 * @param[out] xtangents, ytangents the tangents of the cell faces
 * @param[out] edge_lengths the lengths of the cell faces
 */
typedef void (*fclaw2d_metric_fort_compute_tangents_t)(const int* mx, 
                                                       const int* my, const int* mbc,
                                                       double xd[], double yd[], 
                                                       double zd[],
                                                       double xtangents[], 
                                                       double ytangents[],
                                                       double edge_lengths[]);
    
/**
 * @brief Compute the normals for the cell surfaces
 *
 * @param[in] mx, my the number of cells in the x and y direcitons
 * @param[in] mbc the number of ghost cells
 * @param[in] xlower, ylower the lower lef coordinate of the patch
 * @param[in] xnormals, ynormals the normals of the cell faces
 * @param[in] edge_lengths the lengths of the cell faces
 * @param[in] curvature
 * @param[out] surfnormals the normals of the cell surfaces
 * @param[in] area the area of each cell
 */
typedef void (*fclaw2d_metric_fort_compute_surf_normals_t)(const int* mx, 
                                                           const int* my, const int* mbc,
                                                           double xnormals[], 
                                                           double ynormals[],
                                                           double edge_lengths[],
                                                           double curvature[],
                                                           double surfnormals[], double area[]);


#if 0
typedef double (*fclaw2d_fort_aux_func_t)(double* xc, double *yc);
#endif


/** Fortran subroutine name */
#define FCLAW2D_METRIC_FORT_AVERAGE_AREA  \
         FCLAW_F77_FUNC(fclaw2d_metric_fort_average_area, \
                        FCLAW2D_METRIC_FORT_AVERAGE_AREA) 
/** @copydoc fclaw2d_fort_average_area() */
void FCLAW2D_METRIC_FORT_AVERAGE_AREA(const int* mx, const int* my,
                                 const int* mbc,
                                 double areacoarse[],double areafine[],
                                 const int* igrid);

/** Fortran subroutine name */
#define FCLAW2D_METRIC_FORT_COMPUTE_MESH \
                   FCLAW_F77_FUNC(fclaw2d_metric_fort_compute_mesh, \
                                  FCLAW2D_METRIC_FORT_COMPUTE_MESH)
/** @copydoc fclaw2d_fort_compute_mesh() */
void FCLAW2D_METRIC_FORT_COMPUTE_MESH(const int* mx, const int* my, const int* mbc,
                                      const double* xlower, const double* ylower,
                                      const double* dx, const double* dy,
                                      int* blockno,
                                      double xp[], double yp[], double zp[],
                                      double xd[], double yd[], double zd[]);

/** Fortran subroutine name */
#define FCLAW2D_METRIC_FORT_COMPUTE_AREA \
              FCLAW_F77_FUNC(fclaw2d_metric_fort_compute_area, \
                             FCLAW2D_METRIC_FORT_COMPUTE_AREA)
/** @copydoc fclaw2d_metric_fort_compute_area() */
void FCLAW2D_METRIC_FORT_COMPUTE_AREA(const int* mx, const int* my, const int* mbc,
                               const double* dx, const double* dy,
                               const double* xlower, const double* ylower,
                               const int* blockno, double area[],
                               const int* quadsize, double quadstore[],
                               const int* ghost_only);


/** Fortran subroutine name */
#define FCLAW2D_METRIC_FORT_COMPUTE_NORMALS \
        FCLAW_F77_FUNC(fclaw2d_metric_fort_compute_normals, \
                       FCLAW2D_METRIC_FORT_COMPUTE_NORMALS)

/** @copydoc fclaw2d_metric_fort_compute_normals() */
void FCLAW2D_METRIC_FORT_COMPUTE_NORMALS(const int* mx, const int* my, const int* mbc,
                                  double xp[], double yp[], double zp[],
                                  double xd[], double yd[], double zd[],
                                  double xnormals[], double ynormals[]);

/** Fortran subroutine name */
#define FCLAW2D_METRIC_FORT_COMPUTE_TANGENTS \
            FCLAW_F77_FUNC(fclaw2d_metric_fort_compute_tangents,  \
                           FCLAW2D_METRIC_FORT_COMPUTE_TANGENTS)

/** @copydoc fclaw2d_metric_fort_compute_tangents() */
void FCLAW2D_METRIC_FORT_COMPUTE_TANGENTS(const int* mx, const int* my, const int* mbc,
                                   double xd[], double yd[], double zd[],
                                   double xtangents[], double ytangents[],
                                   double edge_lengths[]);

/** Fortran subroutine name */
#define FCLAW2D_METRIC_FORT_COMPUTE_SURF_NORMALS \
                 FCLAW_F77_FUNC(fclaw2d_metric_fort_compute_surf_normals,  \
                                FCLAW2D_METRIC_FORT_COMPUTE_SURF_NORMALS)

/** @copydoc fclaw2d_metric_fort_compute_surf_normals() */
void FCLAW2D_METRIC_FORT_COMPUTE_SURF_NORMALS(const int* mx, const int* my, 
                                              const int* mbc,
                                              double xnormals[],double ynormals[],
                                              double edge_lengths[],
                                              double curvature[],
                                              double surfnormals[], double area[]);

#ifdef __cplusplus
}
#endif

#endif
