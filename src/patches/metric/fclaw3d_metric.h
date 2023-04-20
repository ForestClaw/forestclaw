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
/**
 * @file 
 * @brief Patch subroutines for patch-related metric terms
 */


#ifndef FCLAW3D_METRIC_H
#define FCLAW3D_METRIC_H

#include <fclaw2d_patch.h>                /* Needed to get enum for build modes */

#include "fclaw3d_metric_default_fort.h"  /* Needed for fort typdefs in vtable */


#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


/** Typedef for ::fclaw3d_metric_vtable */
typedef struct fclaw3d_metric_vtable fclaw3d_metric_vtable_t;

struct fclaw2d_global;
struct fclaw2d_patch;

/* --------------------------- Metric routines (typedefs) ----------------------------- */
/**
 * @brief Compute the cell center and node coordinates for a patch
 *
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw3d_metric_compute_mesh_t)(struct fclaw2d_global *glob,
                                              struct fclaw2d_patch *patch,
                                              int blockno,
                                              int patchno);

/**
 * @brief Compute the area for each cell
 *
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw3d_metric_compute_volume_t)(struct fclaw2d_global *glob,
                                                struct fclaw2d_patch *patch,
                                                int blockno,
                                                int patchno);

/**
 * @brief Compute the volume for each ghost cell
 *
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw3d_metric_compute_volume_ghost_t)(struct fclaw2d_global *glob,
                                                      struct fclaw2d_patch *patch,
                                                      int blockno,
                                                      int patchno);

/**
 * @brief Compute the rotation matrix
 *
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */

typedef void (*fclaw3d_metric_compute_basis_t)(struct fclaw2d_global *glob,
                                               struct fclaw2d_patch *patch,
                                               int blockno,
                                               int patchno);


/* ------------------------------ Metric routines ------------------------------------- */

/**
 * @brief Sets the grid data for a patch and initialized storage.
 * 
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] mx, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] dx, dy the spacings in the x and y directions
 * @param[in] xlower, ylower the lower left coordinate
 * @param[in] xupper, yupper the upper right coordinate
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in] build_mode the build mode
 */

void fclaw3d_metric_patch_define(struct fclaw2d_global* glob,
                                  struct fclaw2d_patch *patch,
                                  int mx, int my, int mz, 
                                  int mbc, 
                                  double dx, double dy, double dz,
                                  double xlower, double ylower, double zlower,
                                  double xupper, double yupper, double zupper,
                                  int blockno, int patchno,
                                  fclaw2d_build_mode_t build_mode);


/**
 * @brief Computes the mesh and tensor arrays
 * 
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 */
void fclaw3d_metric_patch_build(struct fclaw2d_global* glob,
								struct fclaw2d_patch* patch,
								int blockno,
								int patchno);

/**
 * @brief Average the area for coarse patch ghost cells from the fine patches
 * 
 * @param[in] glob the global context
 * @param[in] fine_patches the quad of patchs
 * @param[in,out] coarse_coarse the coarse patch
 * @param[in] blockno the block number
 * @param[in] coarse_patchno the block number of the coarse patch
 * @param[in] fine0_patchno the patch number of the first fine patch
 */
void fclaw3d_metric_patch_build_from_fine(struct fclaw2d_global *glob,
										  struct fclaw2d_patch *fine_patches,
										  struct fclaw2d_patch *coarse_coarse,
										  int blockno,
										  int coarse_patchno,
										  int fine0_patchno);
/**
 * @brief Computes the area array
 * 
 * @param[in] glob the global context
 * @param[in,out] patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 */
void fclaw3d_metric_patch_compute_volume(struct fclaw2d_global *glob,
                                         struct fclaw2d_patch* patch,
                                         int blockno, int patchno);

/* --------------------------------- Access functions --------------------------------- */

/**
 * @brief Get the grid related data for a patch.
 * 
 * @param[in] glob the global context
 * @param[in] patch the patch context
 * @param[out] mx, my the number of cells in the x and y directions
 * @param[out] mbc the number of ghost cells
 * @param[out] xlower, ylower the coordinate of the lower left corner
 * @param[out] dx, dy the spacings in the x and y directions
 */
void fclaw3d_metric_patch_grid_data(struct fclaw2d_global* glob,
									struct fclaw2d_patch* patch,
									int* mx, int* my, int* mz, 
                                    int* mbc,
									double* xlower, double* ylower, double* zlower,
									double* dx, double* dy, double* dz);

/**
 * @brief Get the scalar metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] area the area of each cell in the patch
 * @param[out] face areas the three faces in each cell. 
 *             An array of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:my+mbc+2,3). 
 *             (i,j,k,1) contains the area of face with normal (1,0,0) (left face)
 *             (i,j,k,2) contains the area of face with normal (0,1,0) (front face)
 *             (i,j,k,3) contains the area of face with normal (0,0,1) (bottom face)
 */
void fclaw3d_metric_patch_scalar(struct fclaw2d_global* glob,
								 struct fclaw2d_patch* patch,
								 double **volume, double** faceareas);

/**
 * @brief Get the rotation metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[in]  rotation matrix for face (1,2,3) at cell (i,j,k)
 *             Arrays of dimensions of (-mbc:mx+mbc+2,-mbc:my+mbc+2,N)
 *             where 'N' is user-defined value for the number of components
 *             to store (maximum is 3x(3x3) = 27)
 */
void fclaw3d_metric_patch_basis(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* patch,
                                 double **xrot, double **yrot, double **zrot);

/**
 * @brief Get the mesh metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xp, yp, zp the coordinates of the cell centers
 * @param[out] xd, yd, zd the coordinates of the nodes
 * @param[out] area the area of each cell
 */
void fclaw3d_metric_patch_mesh_data(struct fclaw2d_global* glob,
									struct fclaw2d_patch* patch,
									double **xp, double **yp, double **zp,
									double **xd, double **yd, double **zd,
									double **volume, double** faceareas);

#if 0
/**
 * @brief Get the mesh metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[in]  face areas for each cell
 *             An array of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,3). 
 *             (i,j,k,1) contains the area of face with normal (1,0,0) (left face)
 *             (i,j,k,2) contains the area of face with normal (0,1,0) (front face)
 *             (i,j,k,3) contains the area of face with normal (0,0,1) (bottom face)
 */
void fclaw3d_metric_patch_mesh_data2(struct fclaw2d_global* glob,
									 struct fclaw2d_patch* patch,
                                     double **xrot, double **yrot, double **zrot);
#endif

/**
 * @brief Get the area for each cell of the patch
 * 
 * @param global the global context
 * @param patch the patch to get the are for
 * @return double* the area array
 */
double* fclaw3d_metric_patch_get_volume(struct fclaw2d_global* glob,
                                        struct fclaw2d_patch* patch);


/* ---------------------------- Metric default (virtualized) -------------------------- */

/**
 * @brief @copybrief ::fclaw3d_metric_compute_area_t
 * 
 * Default implementation. Calls fclaw2d_fort_compute_area()
 * 
 * @details @copydetails ::fclaw3d_metric_compute_area_t
 */
void fclaw3d_metric_compute_volume_default(struct fclaw2d_global *glob,
                                           struct fclaw2d_patch* patch,
                                           int blockno, int patchno);


/**
 * @brief @copybrief ::fclaw3d_metric_compute_area_ghost_t
 * 
 * Default implementation. Calls fclaw2d_fort_compute_area()
 * 
 * @details @copydetails ::fclaw3d_metric_compute_area_ghost_t
 */
void fclaw3d_metric_compute_volume_ghost_default(struct fclaw2d_global* glob,
                                                 struct fclaw2d_patch* patch,
                                                 int blockno,
                                                 int patchno);

/**
 * @brief @copybrief ::fclaw3d_metric_compute_mesh_t
 * 
 * Default implementation. Calls fclaw3d_metric_vtable.fort_compute_normals, 
 *                               fclaw3d_metric_vtable.fort_compute_surf_normals, 
 *                               and fclaw3d_metric_vtable.fort_compute_tangents.
 * 
 * @details @copydetails ::fclaw3d_metric_compute_mesh_t
 */
void fclaw3d_metric_compute_mesh_default(struct fclaw2d_global *glob,
										 struct fclaw2d_patch *patch,
										 int blockno,
										 int patchno);


/**
 * @brief @copybrief ::fclaw3d_metric_compute_tensors_t
 * 
 * Default implementation. Calls fclaw3d_metric_vtable.fort_compute_mesh
 * 
 * @details @copydetails ::fclaw3d_metric_compute_tensors_t
 */

void fclaw3d_metric_compute_basis_default(struct fclaw2d_global *glob,
                                          struct fclaw2d_patch *patch,
                                          int blockno,
                                          int patchno);

/* ------------------------------------ Virtual table --------------------------------- */

/**
 * @brief vtable for metric terms
 */
struct fclaw3d_metric_vtable
{
	/* Building patches, including functions to create metric terms */
	/** Computes the mesh coordinates. By default, calls fort_compute_mesh */
	fclaw3d_metric_compute_mesh_t        compute_mesh;

	/** Computes the volume and face area of each cell */
	fclaw3d_metric_compute_volume_t        compute_volume;  /* wrapper */

    /** Computes the area of each ghost cell */
    fclaw3d_metric_compute_volume_ghost_t  compute_volume_ghost;

    /** Computes the basis of each cell */
    fclaw3d_metric_compute_basis_t        compute_basis;  /* wrapper */


	/* ------------------------- Fortran files -----------------------------------------*/
	/** Compute the mesh coordinates */
	fclaw3d_fort_compute_mesh_t         fort_compute_mesh;

    /** Compute the face volume and face areas */
    fclaw3d_fort_compute_volume_t   fort_compute_volume;

    /** Compute the face normal plus two additional vectors in face plane 
        to form a local, orthogonal coordinate system.  Rotation matrix data
        is stored at each of the three faces.  */
    fclaw3d_fort_compute_basis_t       fort_compute_basis;


	/** True if vtable has been set */
	int is_set;
};

/**
 * @brief Get the global vtable variable
 * 
 * @return fclaw3d_metric_vtable_t* the vtable
 */
fclaw3d_metric_vtable_t* fclaw3d_metric_vt(struct fclaw2d_global* glob);

/**
 * @brief Initializes a global vtable variable
 */
void fclaw3d_metric_vtable_initialize(struct fclaw2d_global* glob);


int fclaw3d_metric_patch_nodes_size(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* patch);


#ifdef __cplusplus
}
#endif


#endif /* !FCLAW3D_METRIC_H */
