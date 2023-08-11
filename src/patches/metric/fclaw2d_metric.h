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


#ifndef FCLAW2D_METRIC_H
#define FCLAW2D_METRIC_H

#include <fclaw2d_patch.h>                /* Needed to get enum for build modes */

#include "fclaw2d_metric_default_fort.h"  /* Needed for fort typdefs in vtable */


#ifdef __cplusplus
extern "C"
{
#endif

/** Typedef for ::fclaw2d_metric_vtable */
typedef struct fclaw2d_metric_vtable fclaw2d_metric_vtable_t;

struct fclaw_global;
struct fclaw_patch;

/* --------------------------- Metric routines (typedefs) ----------------------------- */
/**
 * @brief Compute the cell center and node coordinates for a patch
 *
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw2d_metric_compute_mesh_t)(struct fclaw_global *glob,
											  struct fclaw_patch *this_patch,
											  int blockno,
											  int patchno);

/**
 * @brief Compute the area for each cell
 *
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw2d_metric_compute_area_t)(struct fclaw_global *glob,
											  struct fclaw_patch *patch,
											  int blockno,
											  int patchno);

/**
 * @brief Compute the area for each ghost cell
 *
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw2d_metric_compute_area_ghost_t)(struct fclaw_global *glob,
													struct fclaw_patch *patch,
													int blockno,
													int patchno);

/**
 * @brief Compute the basis (normals and tangents)
 *
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] block_no the block number
 * @param[in] patchno the patch number
 */
typedef void (*fclaw2d_metric_compute_basis_t)(struct fclaw_global *glob,
												 struct fclaw_patch *this_patch,
												 int blockno,
												 int patchno);


/* ------------------------------ Metric routines ------------------------------------- */

/**
 * @brief Sets the grid data for a patch and initialized storage.
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] mx, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] dx, dy the spacings in the x and y directions
 * @param[in] xlower, ylower the lower left coordinate
 * @param[in] xupper, yupper the upper right coordinate
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[in] build_mode the build mode
 */
void fclaw2d_metric_patch_define(struct fclaw_global* glob,
								 struct fclaw_patch *this_patch,
								 int mx, int my, int mbc, 
								 double dx, double dy, 
								 double xlower, double ylower,
								 double xupper, double yupper,
								 int blockno, int patchno,
								 fclaw2d_build_mode_t build_mode);


/**
 * @brief Computes the mesh and tensor arrays
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 */
void fclaw2d_metric_patch_build(struct fclaw_global* glob,
								struct fclaw_patch* this_patch,
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
void fclaw2d_metric_patch_build_from_fine(struct fclaw_global *glob,
										  struct fclaw_patch *fine_patches,
										  struct fclaw_patch *coarse_coarse,
										  int blockno,
										  int coarse_patchno,
										  int fine0_patchno);
/**
 * @brief Computes the area array
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 */
void fclaw2d_metric_patch_compute_area(struct fclaw_global *glob,
									   struct fclaw_patch* this_patch,
									   int blockno, int patchno);

/* --------------------------------- Access functions --------------------------------- */

/**
 * @brief Get the grid related data for a patch.
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[out] mx, my the number of cells in the x and y directions
 * @param[out] mbc the number of ghost cells
 * @param[out] xlower, ylower the coordinate of the lower left corner
 * @param[out] dx, dy the spacings in the x and y directions
 */
void fclaw2d_metric_patch_grid_data(struct fclaw_global* glob,
									struct fclaw_patch* this_patch,
									int* mx, int* my, int* mbc,
									double* xlower, double* ylower,
									double* dx, double* dy);

/**
 * @brief Get the scalar metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] area the area of each cell in the patch
 * @param[out] edgelengths the edge lenghts for each cell. 
 *             An array of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the the bottom edge length for cell (i,j).
 *             (i,j,2) contains the the left edge length for cell (i,j).
 * @param[out] curvature the curvature for each cell in the patch
 */
void fclaw2d_metric_patch_scalar(struct fclaw_global* glob,
								 struct fclaw_patch* this_patch,
								 double **area, double** edgelengths,
								 double **curvature);

/**
 * @brief Get the vector metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xnormals, ynormals the normals of each face of the patch
 *             Arrays of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the normal for the bottom edge of cell (i,j).
 *             (i,j,2) contains the normal for the left edge of cell (i,j).
 * @param[out] xtangents, ytangents the tangents of each face of the patch
 *             Arrays of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the tangent for the bottom edge of cell (i,j).
 *             (i,j,2) contains the tangent for the left edge of cell (i,j).
 * @param[out] surfnormals the surface normal for each cell center of the patch
 *             An array of dimension(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)
 */
void fclaw2d_metric_patch_vector(struct fclaw_global* glob,
                                 struct fclaw_patch* this_patch,
                                 double **xnormals, double **ynormals,
                                 double **xtangents, double **ytangents,
                                 double **surfnormals);

/**
 * @brief Get the mesh metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xp, yp, zp the coordinates of the cell centers
 * @param[out] xd, yd, zd the coordinates of the nodes
 * @param[out] area the area of each cell
 */
void fclaw2d_metric_patch_mesh_data(struct fclaw_global* glob,
									struct fclaw_patch* this_patch,
									double **xp, double **yp, double **zp,
									double **xd, double **yd, double **zd,
									double **area);

/**
 * @brief Get the mesh metrics of a patch
 * 
 * @param[in]  glob glob the global context
 * @param[in]  patch the patch context
 * @param[out] xnormals, ynormals the normals of each face of the patch
 *             Arrays of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the normal for the bottom edge of cell (i,j).
 *             (i,j,2) contains the normal for the left edge of cell (i,j).
 * @param[out] xtangents, ytangents the tangents of each face of the patch
 *             Arrays of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the tangent for the bottom edge of cell (i,j).
 *             (i,j,2) contains the tangent for the left edge of cell (i,j).
 * @param[out] surfnormals the surface normal for each cell center of the patch
 *             An array of dimension(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)
 * @param[out] edgelengths the edge lenghts for each cell. 
 *             An array of dimension (-mbc:mx+mbc+2,-mbc:my+mbc+2,2). 
 *             (i,j,1) contains the the bottom edge length for cell (i,j).
 *             (i,j,2) contains the the left edge length for cell (i,j).
 * @param[out] curvature the curvature for each cell in the patch
 */
void fclaw2d_metric_patch_mesh_data2(struct fclaw_global* glob,
									 struct fclaw_patch* this_patch,
									 double **xnormals, double **ynormals,
									 double **xtangents, double **ytangents,
									 double **surfnormals,
									 double **edgelengths, double **curvature);

/**
 * @brief Get the area for each cell of the patch
 * 
 * @param global the global context
 * @param this_patch the patch to get the are for
 * @return double* the area array
 */
double* fclaw2d_metric_patch_get_area(struct fclaw_global* glob,
	                                  struct fclaw_patch* this_patch);


/* ---------------------------- Metric default (virtualized) -------------------------- */

/**
 * @brief @copybrief ::fclaw2d_metric_compute_area_t
 * 
 * Default implementation. Calls fclaw2d_fort_compute_area()
 * 
 * @details @copydetails ::fclaw2d_metric_compute_area_t
 */
void fclaw2d_metric_compute_area_default(struct fclaw_global *glob,
										 struct fclaw_patch* this_patch,
										 int blockno, int patchno);


/**
 * @brief @copybrief ::fclaw2d_metric_compute_area_ghost_t
 * 
 * Default implementation. Calls fclaw2d_fort_compute_area()
 * 
 * @details @copydetails ::fclaw2d_metric_compute_area_ghost_t
 */
void fclaw2d_metric_compute_area_ghost_default(struct fclaw_global* glob,
											   struct fclaw_patch* this_patch,
											   int blockno,
											   int patchno);

/**
 * @brief @copybrief ::fclaw2d_metric_compute_mesh_t
 * 
 * Default implementation. Calls fclaw2d_metric_vtable.fort_compute_normals, 
 *                               fclaw2d_metric_vtable.fort_compute_surf_normals, 
 *                               and fclaw2d_metric_vtable.fort_compute_tangents.
 * 
 * @details @copydetails ::fclaw2d_metric_compute_mesh_t
 */
void fclaw2d_metric_compute_mesh_default(struct fclaw_global *glob,
										 struct fclaw_patch *this_patch,
										 int blockno,
										 int patchno);


/**
 * @brief @copybrief ::fclaw2d_metric_compute_basis_t
 * 
 * Default implementation. Calls fclaw2d_metric_vtable.fort_compute_mesh
 * 
 * @details @copydetails ::fclaw2d_metric_compute_basis_t
 */
void fclaw2d_metric_compute_basis_default(struct fclaw_global *glob,
											struct fclaw_patch *this_patch,
											int blockno,
											int patchno);


/* ------------------------------------ Virtual table --------------------------------- */

/**
 * @brief vtable for metric terms
 */
struct fclaw2d_metric_vtable
{
	/* Building patches, including functions to create metric terms */
	/** Computes the mesh coordinates. By default, calls fort_compute_mesh */
	fclaw2d_metric_compute_mesh_t        compute_mesh;
	/** Computes the area of each cell */
	fclaw2d_metric_compute_area_t        compute_area;  /* wrapper */
	/** Computes the area of each ghost cell */
	fclaw2d_metric_compute_area_ghost_t  compute_area_ghost;
	/** Computes the basis. By default, calls fort_compute_normals, fort_compute_tangents, and fort_compute_surf_normals */
	fclaw2d_metric_compute_basis_t     compute_basis;  /* wrapper */

	/* Fortran files */
	/** Compute the mesh coordinates */
	fclaw2d_metric_fort_compute_mesh_t          fort_compute_mesh;
	/** Compute the face normals */
	fclaw2d_metric_fort_compute_normals_t       fort_compute_normals;
	/** Compute the face tangents */
	fclaw2d_metric_fort_compute_tangents_t      fort_compute_tangents;
	/** Compute the surface normals */
	fclaw2d_metric_fort_compute_surf_normals_t  fort_compute_surf_normals;

	/** True if vtable has been set */
	int is_set;
};

/**
 * @brief Get the global vtable variable
 * 
 * @return fclaw2d_metric_vtable_t* the vtable
 */
fclaw2d_metric_vtable_t* fclaw2d_metric_vt(struct fclaw_global* glob);

/**
 * @brief Initializes a global vtable variable
 */
void fclaw2d_metric_vtable_initialize(struct fclaw_global* glob);

int fclaw2d_metric_patch_nodes_size(struct fclaw_global* glob,
                                    struct fclaw_patch* patch);



#ifdef __cplusplus
}
#endif


#endif /* !FCLAW2D_METRIC_H */
