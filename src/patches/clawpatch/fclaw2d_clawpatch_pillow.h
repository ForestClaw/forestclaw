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

#ifndef PILLOWSPHERE_H
#define PILLOWSPHERE_H

#include <fclaw2d_defs.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_glob;
struct fclaw2d_patch;
struct fclaw2d_patch_transform;

/**
 * @file
 * @brief This handles the boundary conditions at the block
 * corners for the pillow sphere.
 */

/** typedef */
typedef struct fclaw2d_clawpatch_pillow_vtable fclaw2d_clawpatch_pillow_vtable_t;


/* ----------------------------- Fortran typedefs ------------------------------------- */

/**
 * @brief Handles the boundary condition at block corners
 * 
 * @param[in] my, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] meqn the number of equations
 * @param[in,out] qthis this solution
 * @param[in,out] qneighbor the neighbor solution
 * @param[in] icorner the corner that the neighbor is on
 * @param[in] iblock the block number
 */
typedef void (*pillow_fort_copy_block_corner_t)(int* mx, int* my,
                                                int* mbc, int* meqn,
                                                double qthis[], 
                                                double qneighbor[], 
                                                int* icorner,
                                                int* iblock);

/**
 * @brief Handles the boundary condition at block corners with finer neighbors
 * 
 * @param[in] my, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] meqn the number of equations
 * @param[in] refratio the refinement ratio
 * @param[in,out] qcoarse this solution
 * @param[in] qfine the neighbor solution
 * @param[in] areacoarse cell areas
 * @param[in] areafine neighbor cell areas
 * @param[in] icorner the corner that the fine neighbor is on
 * @param[in] blockno the block number
 */
typedef void  (*pillow_fort_average_block_corner_t)(int* mx, int* my, int* mbc,
                                                    int* meqn, 
                                                    int* refratio, 
                                                    double qcoarse[],
                                                    double qfine[], 
                                                    double areacoarse[], 
                                                    double areafine[],
                                                    int* icorner,
                                                    int* blockno);

/**
 * @brief Handles the boundary condition at block corners with coarser neighbors
 * 
 * @param[in] my, my the number of cells in the x and y directions
 * @param[in] mbc the number of ghost cells
 * @param[in] meqn the number of equations
 * @param[in] refratio the refinement ratio
 * @param[in] qcoarse this solution
 * @param[in,out] qfine the neighbor solution
 * @param[in] icorner_coarse the corner that the fine neighbor is on
 * @param[in] blockno the block number
 */
typedef void  (*pillow_fort_interpolate_block_corner_t)(int* mx, int* my, int* mbc,
                                                        int* meqn, int* refratio,
                                                        double qcoarse[],
                                                        double qfine[], 
                                                        int* icoarse_corner,
                                                        int* blockno);
#if 0
typedef void (*pillow_fort_copy_block_corner_t)(const int* mx, 
                                                const int* my,
                                                const int *mz,
                                                const int* mbc, 
                                                const int* meqn,
                                                double qthis[], 
                                                double qneighbor[], 
                                                const int* icorner,
                                                const int* iblock);

typedef void  (*pillow_fort_average_block_corner_t)(const int* mx, 
                                                    const int* my, 
                                                    const int* mz, 
                                                    const double *dz, 
                                                    const int* mbc,
                                                    const int* meqn, 
                                                    const int* refratio, 
                                                    double qcoarse[],
                                                    double qfine[], 
                                                    double areacoarse[], 
                                                    double areafine[],
                                                    const int* coarse_corner,
                                                    const int* blockno);

typedef void  (*pillow_fort_interpolate_block_corner_t)(const int* mx, 
                                                        const int* my, 
                                                        const int* mz,
                                                        const int* mbc,
                                                        const int* meqn, 
                                                        const int* refratio,
                                                        double qcoarse[],
                                                        double qfine[], 
                                                        const int* icoarse_corner,
                                                        const int* blockno);

#endif

/* ----------------------------- Use pillow sphere ------------------------------------ */

/**
 * @brief Sets global patch_vtable to use pollow sphere routines
 */
void fclaw2d_clawpatch_use_pillowsphere();

/* --------------------------------- Virtual table ------------------------------------ */

/**
 * @brief Initialize the pillow vtable
 * 
 * @param glob the global context
 * @param claw_version the clawaptck verstion 4 for v4.6, 5 for v5
 */
void fclaw2d_clawpatch_pillow_vtable_initialize(struct fclaw2d_global* glob, 
                                                int claw_version);

/**
 * @brief vtable for handling block corners for pillow sphere
 */
struct fclaw2d_clawpatch_pillow_vtable
{
    /* Block corners */
    /** Handles the boundary condition at block corners */
    pillow_fort_copy_block_corner_t           fort_copy_block_corner;
    /** Handles the boundary condition at block corners with finer neighbors */
    pillow_fort_average_block_corner_t        fort_average_block_corner;
    /** Handles the boundary condition at block corners with coarser neighbors */
    pillow_fort_interpolate_block_corner_t    fort_interpolate_block_corner;

    /** True if vtable is set */
    int is_set;
};


/* ------------------------------- Public access functions ---------------------------- */

/**
 * @brief Get the pillow vtable
 * 
 * @param glob the global context
 * @return fclaw2d_clawpatch_pillow_vtable_t* the vtable
 */
fclaw2d_clawpatch_pillow_vtable_t* fclaw2d_clawpatch_pillow_vt(struct fclaw2d_global* glob);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif    