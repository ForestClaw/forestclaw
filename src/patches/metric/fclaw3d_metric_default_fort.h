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

#ifndef FCLAW3D_METRIC_DEFAULT_FORT_H
#define FCLAW3D_METRIC_DEFAULT_FORT_H

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
typedef void (*fclaw3d_fort_compute_mesh_t)(const int* mx, const int* my, const int *mz,
                                            const int* mbc,
                                            const double* xlower, const double* ylower,
                                            const double *zlower,
                                            const double* dx, const double* dy,
                                            const double *dz,
                                            int* blockno,
                                            double xp[], double yp[], double zp[],
                                            double xd[], double yd[], double zd[]);

/* This doesn't appear to be used ... */
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
typedef void (*fclaw3d_fort_compute_volume_t)(const int* mx, const int* my, const int* mz,
                                              const int* mbc,
                                              const double* dx, const double* dy,
                                              const double *dz,
                                              const double* xlower, const double* ylower,
                                              const double* zlower,
                                              const int* blockno, 
                                              double xd[], double yd[], double zd[],
                                              double volume[], double faceareas[],
                                              const int* hexsize, double* hexfine,
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
typedef void (*fclaw3d_fort_compute_basis_t)(const int* mx, const int* my, const int* mz, 
                                             const int* mbc,
                                             double xd[], double yd[], double zd[],
                                             double xrot[], double yrot[], double zrot[],
                                             const int* ghost_only);

/** Fortran subroutine name */
#define FCLAW3D_METRIC_FORT_COMPUTE_MESH \
                  FCLAW_F77_FUNC(fclaw3d_metric_fort_compute_mesh, \
                                 FCLAW3D_METRIC_FORT_COMPUTE_MESH)
/** @copydoc fclaw3d_fort_compute_mesh() */
void FCLAW3D_METRIC_FORT_COMPUTE_MESH(const int* mx, const int* my, const int* mz, 
                               const int* mbc,
                               const double* xlower, const double* ylower,
                               const double* zlower, 
                               const double* dx, const double* dy, 
                               const double *dz,
                               int* blockno,
                               double xp[], double yp[], double zp[],
                               double xd[], double yd[], double zd[]);

/** Fortran subroutine name */
#define FCLAW3D_METRIC_FORT_COMPUTE_VOLUME \
                FCLAW_F77_FUNC(fclaw3d_metric_fort_compute_volume, \
                               FCLAW3D_METRIC_FORT_COMPUTE_VOLUME)
/** @copydoc fclaw3d_fort_compute_volume() */
void FCLAW3D_METRIC_FORT_COMPUTE_VOLUME(const int* mx, const int* my, const int* mz,
                                 const int* mbc,
                                 const double* dx, const double* dy,
                                 const double* dz,
                                 const double* xlower, const double* ylower,
                                 double *zlower,
                                 const int* blockno, 
                                 double xd[], double yd[], double zd[],
                                 double volume[],
                                 double faceareas[],
                                 const int* hexsize, double hexfine[],
                                 const int* ghost_only);



/** Fortran subroutine name */
#define FCLAW3D_METRIC_FORT_COMPUTE_BASIS \
                       FCLAW_F77_FUNC(fclaw3d_metric_fort_compute_basis, \
                                      FCLAW3D_METRIC_FORT_COMPUTE_BASIS)

/** @copydoc fclaw3d_fort_compute_tangents() */
void FCLAW3D_METRIC_FORT_COMPUTE_BASIS(const int* mx, const int* my, const int* mz,
                                const int* mbc,
                                double xd[], double yd[], double zd[],
                                double xrot[], double yrot[], double zrot[],
                                const int* ghost_only);

/** Fortran subroutine name */
#define FCLAW3DX_METRIC_FORT_AVERAGE_VOLUME \
                       FCLAW_F77_FUNC(fclaw3dx_metric_fort_average_volume, \
                                      FCLAW3DX_METRIC_FORT_AVERAGE_VOLUME) 
/** @copydoc fclaw3d_fort_average_volume() */
void FCLAW3DX_METRIC_FORT_AVERAGE_VOLUME(const int* mx, const int* my, const int* mz,
                                  const int* mbc,
                                  double volcoarse[],double volfine[],
                                  const int* igrid);

/** Fortran subroutine name */
#define FCLAW3DX_METRIC_FORT_AVERAGE_FACEAREA \
                FCLAW_F77_FUNC(fclaw3dx_metric_fort_average_facearea, \
                               FCLAW3DX_METRIC_FORT_AVERAGE_FACEAREA) 
/** @copydoc fclaw3d_fort_average_volume() */
void FCLAW3DX_METRIC_FORT_AVERAGE_FACEAREA(const int* mx, const int* my, const int* mz,
                                    const int* mbc,
                                    double fa_coarse[],double fa_fine[],
                                    const int* igrid);


#ifdef __cplusplus
}
#endif

#endif
