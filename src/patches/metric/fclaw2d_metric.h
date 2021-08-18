/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW2D_METRIC_H
#define FCLAW2D_METRIC_H

#include <fclaw2d_patch.h>                /* Needed to get enum for build modes */

#include "fclaw2d_metric_default_fort.h"  /* Needed for fort typdefs in vtable */


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef struct fclaw2d_metric_vtable fclaw2d_metric_vtable_t;

struct fclaw2d_global;
struct fclaw2d_patch;

/* --------------------------- Metric routines (typedefs) ----------------------------- */

typedef void (*fclaw2d_metric_compute_mesh_t)(struct fclaw2d_global *glob,
											  struct fclaw2d_patch *this_patch,
											  int blockno,
											  int patchno);

typedef void (*fclaw2d_metric_compute_area_t)(struct fclaw2d_global *glob,
											  struct fclaw2d_patch *patch,
											  int blockno,
											  int patchno);

typedef void (*fclaw2d_metric_compute_area_ghost_t)(struct fclaw2d_global *glob,
													struct fclaw2d_patch *patch,
													int blockno,
													int patchno);

typedef void (*fclaw2d_metric_compute_tensors_t)(struct fclaw2d_global *glob,
												 struct fclaw2d_patch *this_patch,
												 int blockno,
												 int patchno);


/* ------------------------------ Metric routines ------------------------------------- */

void fclaw2d_metric_patch_define(struct fclaw2d_global* glob,
								 struct fclaw2d_patch *this_patch,
								 int mx, int my, int mbc, 
								 double dx, double dy, 
								 double xlower, double ylower,
								 double xupper, double yupper,
								 int blockno, int patchno,
								 fclaw2d_build_mode_t build_mode);


void fclaw2d_metric_patch_setup(struct fclaw2d_global* glob,
								struct fclaw2d_patch* this_patch,
								int blockno,
								int patchno);

void fclaw2d_metric_patch_setup_from_fine(struct fclaw2d_global *glob,
										  struct fclaw2d_patch *fine_patches,
										  struct fclaw2d_patch *coarse_coarse,
										  int blockno,
										  int coarse_patchno,
										  int fine0_patchno);

void fclaw2d_metric_patch_compute_area(struct fclaw2d_global *glob,
									   struct fclaw2d_patch* this_patch,
									   int blockno, int patchno);

/* --------------------------------- Access functions --------------------------------- */

void fclaw2d_metric_patch_grid_data(struct fclaw2d_global* glob,
									struct fclaw2d_patch* this_patch,
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
void fclaw2d_metric_patch_scalar(struct fclaw2d_global* glob,
								 struct fclaw2d_patch* this_patch,
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
void fclaw2d_metric_patch_vector(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch,
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
void fclaw2d_metric_patch_mesh_data(struct fclaw2d_global* glob,
									struct fclaw2d_patch* this_patch,
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
void fclaw2d_metric_patch_mesh_data2(struct fclaw2d_global* glob,
									 struct fclaw2d_patch* this_patch,
									 double **xnormals, double **ynormals,
									 double **xtangents, double **ytangents,
									 double **surfnormals,
									 double **edgelengths, double **curvature);

/**
 * @brief Get the area for each cell of the patch
 * 
 * @param this_patch the patch to get the are for
 * @return double* the area array
 */
double* fclaw2d_metric_patch_get_area(struct fclaw2d_patch* this_patch);


/* ---------------------------- Metric default (virtualized) -------------------------- */

void fclaw2d_metric_compute_area_default(struct fclaw2d_global *glob,
										 struct fclaw2d_patch* this_patch,
										 int blockno, int patchno);


void fclaw2d_metric_compute_area_ghost_default(struct fclaw2d_global* glob,
											   struct fclaw2d_patch* this_patch,
											   int blockno,
											   int patchno);

void fclaw2d_metric_compute_mesh_default(struct fclaw2d_global *glob,
										 struct fclaw2d_patch *this_patch,
										 int blockno,
										 int patchno);


void fclaw2d_metric_compute_tensors_default(struct fclaw2d_global *glob,
											struct fclaw2d_patch *this_patch,
											int blockno,
											int patchno);


/* ------------------------------------ Virtual table --------------------------------- */

fclaw2d_metric_vtable_t* fclaw2d_metric_vt();


struct fclaw2d_metric_vtable
{
	/* Building patches, including functions to create metric terms */
	fclaw2d_metric_compute_mesh_t        compute_mesh;    /* wrapper */
	fclaw2d_metric_compute_area_t        compute_area;  /* wrapper */
	fclaw2d_metric_compute_area_ghost_t  compute_area_ghost;
	fclaw2d_metric_compute_tensors_t     compute_tensors;  /* wrapper */

	/* Fortran files */
	fclaw2d_fort_compute_mesh_t          fort_compute_mesh;
	fclaw2d_fort_compute_normals_t       fort_compute_normals;
	fclaw2d_fort_compute_tangents_t      fort_compute_tangents;
	fclaw2d_fort_compute_surf_normals_t  fort_compute_surf_normals;

	int is_set;
};

void fclaw2d_metric_vtable_initialize();


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FCLAW2D_METRIC_H */
