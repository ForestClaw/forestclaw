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

#ifndef FCLAW3DX_CLAWPATCH_HPP
#define FCLAW3DX_CLAWPATCH_HPP

/**
 * @file 
 * C++ structures for clawpatch
 */
#include <fclaw2d_farraybox.hpp>  /* Needed for FArray boxes */

struct fclaw2d_patch;
struct fclaw2d_global;
struct fclaw2d_metric_patch_t;

/**
 * @brief Stores data for each patch
 */
class fclaw3dx_clawpatch_t
{
public :
    /* Solution data */
    int meqn; /**< number of equations */                   
    FArrayBox griddata; /**< the current solution */
    FArrayBox griddata_last; /**< the solution at the last timestep */
    FArrayBox griddata_save; /**< the saved solution */
    FArrayBox griddata_time_interpolated; /**< the time interpolated solution */
    FArrayBox griderror; /**< the error */

    /** Exact solution for diagnostics */
    FArrayBox exactsolution;

    int mfields;  /**< Number of fields in the rhs */
    FArrayBox rhs;  /**< RHS for elliptic problems */

    FArrayBox elliptic_error;  /**< Error for elliptic problems */
    FArrayBox elliptic_soln;  /**< Solution for elliptic problems */

    /** Registers for accumulating mismatches at coarse/fine interfaces */
    struct fclaw3dx_clawpatch_registers *registers;

    /* Grid info */
    int mx; /**< number of cells in the x direction */          
    int my; /**< number of cells in the y direction */  
    int mz; /**< number of cells in the z direction */  
    int mbc; /**< number ghost cells */
    int maux; /**< number aux equations */

    double dx; /**< cell spacing in the x direction */
    double dy; /**< cell spacing in the y direction */
    double dz; /**< cell spacing in the z direction */

    double xlower; /**< x coordinate of the left face of the patch */
    double ylower; /**< y coordinate of the bottom face of the patch */
    double zlower; /**< z coordinate of the lower face of the patch */
    double xupper; /**< x coordinate of the right face of the patch */
    double yupper; /**< y coordinate of the top face of the patch */
    double zupper; /**< z coordinate of the upper face of the patch */

    /** Auxilliary array (used by Clawpack 4.6 and 5.0) */
    FArrayBox aux;
    FArrayBox aux_save;

    /* Mapping and metric info */
    int manifold; /**< true if using manifold */   
    int blockno; /**< the block number of a patch */

    fclaw2d_metric_patch_t *mp; /**< the metric data for a patch */

    /** Extra storage needed by the solver(s) */
    void* solver_data;

    /** User data*/ 
    void* user_data;
};

/**
 * @brief Get the clawpatch structure for a patch
 * 
 * @param this_patch the patch context
 * @return fclaw3dx_clawpatch_t* the clawpatch structure
 */
fclaw3dx_clawpatch_t* 
fclaw3dx_clawpatch_get_clawpatch(struct fclaw2d_patch* this_patch);

/**
 * @brief Get the metrix structure for a patch
 * 
 * @param this_patch the patch context
 * @return fclaw2d_metric_patch_t* the metric structure
 */
fclaw2d_metric_patch_t* 
fclaw3dx_clawpatch_get_metric_patch(struct fclaw2d_patch* this_patch);




#endif /* !FCLAW2D_CLAWPATCH_HPP */
