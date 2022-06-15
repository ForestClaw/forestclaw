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

#ifndef FCLAW3DX_CLAWPATCH_CONSERVATION_H
#define FCLAW3DX_CLAWPATCH_CONSERVATION_H

#include <fclaw_clawpatch_enums.h>

/** 
 *  @file
 *  Clawpatch conservation related structures and routines
 */

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw2d_patch_transform_data;
struct fclaw2d_global;
struct fclaw2d_patch;

/**
 * @brief fclaw3dx_clawpatch_registers type
 */
typedef struct fclaw3dx_clawpatch_registers fclaw3dx_clawpatch_registers_t;

/**
 * @brief Stores flux information for the edges of patch
 */
struct fclaw3dx_clawpatch_registers
{
	/** Stores the value of the flux function at cell centers along the edges */
	double *edge_fluxes[4];

	/* Scaling factors */
	/** Array of edge lengths for each edge of the patch */
	double *edgelengths[4];
	/** Array of areas for each edge of the patch */
	double *area[4];

	/* 1d arrays stored on left/right faces (fp/fm) and top/bottom faces (gp/gm) */

    /** Fluxes along the left face */
	double *fp[2];
    /** Fluxes along the right face */
	double *fm[2];
    /** Fluxes along the bottom face */
	double *gp[2];
    /** Fluxes along the top face */
	double *gm[2];
};

/**
 * @brief Adds fine grid corrections to coarse grid.  
 *
 * @param[in]     mx, my the number of cells in the x and y directions
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of ghost cells
 * @param[in]     idir the direction of the interface 0 for bottom/top 
 *                1 for left/right
 * @param[in]     iface_coarse the interface on the coarse patch
 * @param[in]     coarse_blockno the block number of the coarse patch
 * @param[in]     fine_blockno the block number of the fine patch
 * @param[in]     normal_match true if normals on coarse and fine patches match
 * @param[in]     area0, area1, area2, area3 area of cells along the 
 *                edges
 * @param[in,out] qcoarse the solution on the coarse patch
 * @param[in]     fpthis0, fmthis1, gpthis2, gmthis3 the fluxes along 
 *                the edges of the coarse patch
 * @param[in]     fmfine0, fmfine1, gpfine2, fpfine3 the fluxes along
 *                the dges of the fine patch
 * @param[in]     efthis0, efthis1, efthis2, efthis3 the values of the
 *                flux function along the coarse patch edges
 * @param[in]     eff0, eff1, eff2, eff3 the values of the flux 
 *                function along the fine patch edges
 * @param[in,out] qfine_dummy work vector for determining values in the fine patch
 * @param[in]     transform_cptr the pointer to the fclaw2d_patch_transform_data struct
 */ 
typedef void  (*fclaw3dx_clawpatch_fort_time_sync_f2c_t)(const int* mx,
                        const int* my,
                        const int *mbc,
                        const int *meqn,
                        const int* idir,
                        const int* iface_coarse,
                        const int* coarse_blockno,
                        const int* fine_blockno,
                        const int* normal_mismatch,
                        double area0[], double area1[],
                        double area2[], double area3[],
                        double qcoarse[], 
                        double fmthis0[], 
                        double fpthis1[],
                        double gmthis2[], 
                        double gpthis3[],
                        double fmfine0[], double fpfine1[],
                        double gmfine2[], double gpfine3[],
                        double efthis0[], double efthis1[],
                        double efthis2[], double efthis3[],
                        double eff0[], double eff1[],
                        double eff2[], double eff3[],
                        double qfine_dummy[],
                        struct fclaw2d_patch_transform_data** 
                        transform_cptr);

/**
 * @brief Adds wave corrections at same level interfaces.  
 * 
 * This accounts for metric mismatches that can occur at block boundaries.
 *
 * @param[in]     mx, my the number of cells in the x and y directions
 * @param[in]     mbc the number of ghost cells
 * @param[in]     meqn the number of ghost cells
 * @param[in]     idir the direction of the interface 0 for bottom/top 
 *                1 for left/right
 * @param[in]     iface_this the interface on the this patch
 * @param[in]     this_blockno the block number of the coarse patch
 * @param[in]     nbr_blockno the block number of the fine patch
 * @param[in]     area0, area1, area2, area3 area of cells along the 
 *                edges on this patch
 * @param[in,out] qthis the solution on this patch
 * @param[in]     fpthis0, fmthis1, gpthis2, gmthis3 the fluxes along 
 *                the edges of the this patch
 * @param[in]     fmfine0, fmfine1, gpfine2, fpfine3 the fluxes along
 *                the dges of the neighbor patch
 * @param[in]     efthis0, efthis1, efthis2, efthis3 the values of the
 *                flux function along the this patch edges
 * @param[in]     efnbr0, efnbr1, efnbr2, efnbr3 the values of the flux 
 *                function along the neighbor patch edges
 * @param[in,out] qnbr_dummy work vector for determining values in the fine patch
 * @param[in]     transform_cptr the pointer to the fclaw2d_patch_transform_data struct
 */
typedef void  (*fclaw3dx_clawpatch_fort_time_sync_samesize_t)(const int* mx,
                                                     const int* my,
                                                     const int *mbc,
                                                     const int *meqn,
                                                     const int* idir,
                                                     const int* iface_this,
                                                     const int* this_blockno,
                                                     const int* nbr_blockno,
                                                     double area0[], double area1[],
                                                     double area2[], double area3[],
                                                     double qthis[], 
                                                     double fmthis0[], 
                                                     double fpthis1[],
                                                     double gmthis2[], 
                                                     double gpthis3[],
                                                     double fmnbr0[], double fpnbr1[],
                                                     double gmnbr2[], double gpnbr3[],
                                                     double efthis0[], double efthis1[],
                                                     double efthis2[], double efthis3[],
                                                     double efnbr0[], double efnbr1[],
                                                     double efnbr2[], double efnbr3[],
                                                     double qnbr_dummy[],
                                                     struct fclaw2d_patch_transform_data** 
                                                     transform_cptr);



/**
 * @brief Allocate a new registers struct
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 * @param[out] registers the newly allocated registers struct
 */
void fclaw3dx_clawpatch_time_sync_new(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* this_patch,
                                     int blockno,int patchno,
                                     fclaw3dx_clawpatch_registers_t **registers);

/**
 * @brief Deallocate the registers struct
 * 
 * @param[in,out] registers the registers, set to NULL on return
 */
void fclaw3dx_clawpatch_time_sync_delete(fclaw3dx_clawpatch_registers_t **registers);

/**
 * @brief Inialize the area and edgelength arrays of the registers
 * 
 * @param[in] glob the global context
 * @param[in] this_patch the patch context
 * @param[in] blockno the block number
 * @param[in] patchno the patch number
 */
void fclaw3dx_clawpatch_time_sync_setup(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch* this_patch,
                                       int blockno,int patchno);

/**
 * @brief Adds fine grid corrections to coarse grid.  
 * 
 * @param[in] glob the global context
 * @param[in,out] coarse_patch the coarse patch context
 * @param[in] fine_patch the fine patch context
 * @param[in] coarse_blockno the block number of the coarse patch
 * @param[in] fine_blockno the block number of the fine patch
 * @param[in] coarse_patchno the patch number of the coarse patch
 * @param[in] idir the direction of the interface 0 for bottom/top 
 *            1 for left/right
 * @param[in] igrid the index of the fine grid in the child array
 * @param[in] iface_coarse the interface on the coarse patch
 * @param[in] time_interp NOT USED
 * @param[in] transform_data the transform for the neighbor's coordinates
 */
void fclaw3dx_clawpatch_time_sync_f2c(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* coarse_patch,
                                     struct fclaw2d_patch* fine_patch,
                                     int coarse_blockno, int fine_blockno,
                                     int coarse_patchno, 
                                     int idir,
                                     int igrid,
                                     int iface_coarse,
                                     int time_interp,
                                     struct fclaw2d_patch_transform_data
                                     *transform_data);
	
/**
 * @brief Adds corrections to patches that are at the same levle and are at block boundaries.
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch this patch
 * @param[in] neighbor_patch the neighbor patch
 * @param[in] this_iface the interface that the neighbor patch is on
 * @param[in] idir the direction of the interface 0 for bottom/top 
 *            1 for left/right
 * @param[in] transform_data the transform for the neighbor's coordinates
 */
void fclaw3dx_clawpatch_time_sync_samesize(struct fclaw2d_global* glob,
                                          struct fclaw2d_patch* this_patch,
                                          struct fclaw2d_patch* neighbor_patch,
                                          int this_iface,int idir,
                                          struct fclaw2d_patch_transform_data
                                          *transform_data);

/**
 * @brief Resets arrays in the registers
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in] coarse_level the the level of the coarse patch
 * @param[in] reset_mode the reset mode ::fclaw2d_time_sync_type
 */
void fclaw3dx_clawpatch_time_sync_reset(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch *this_patch,
                                       int coarse_level,
                                       int reset_mode);

/**
 * @brief Packs/Unpacks the registers to/from a buffer
 * 
 * @param[in] glob the global context
 * @param[in,out] this_patch the patch context
 * @param[in,out] qpack the buffer
 * @param[in] frsize the size of the buffer
 * @param[in] packmode the packing mode
 * @param[out] ierror the error value
 */
void fclaw3dx_clawpatch_time_sync_pack_registers(struct fclaw2d_global *glob,
                                                struct fclaw2d_patch *this_patch,
                                                double *qpack,
                                                int frsize, 
                                                fclaw_clawpatch_packmode_t packmode, 
                                                int *ierror);

#ifdef __cplusplus
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_CONSERVATION_H */

