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
 * C++ structures for metric patch
 */

#ifndef FCLAW2D_METRIC_HPP
#define FCLAW2D_METRIC_HPP

#include <fclaw_farraybox.hpp>     /* Needed for FArray box used to store metric data */

/**
 * @brief Struct for patch metric data
 */
struct fclaw2d_metric_patch_t
{
    /** The number of cells in the x direction */
    int mx;           
    /** The number of cells in the x direction */
    int my;           
    /** The number of ghost cells */
    int mbc;
    /** The block number */
    int blockno;

    /** The spacing in the x direction */
    double dx;
    /** The spacing in the y direction */
    double dy;
    /** x value of bottom face */
    double xlower;
    /** y value of left face */
    double ylower;
    /** x value of upper face */
    double xupper;
    /** y value of right face */
    double yupper;
    
    /** x coordinates of cell centers */
    FArrayBox xp;
    /** y coordinates of cell centers */
    FArrayBox yp;
    /** z coordinates of cell centers */
    FArrayBox zp;

    /** x coordinates of nodes */
    FArrayBox xd;
    /** y coordinates of nodes */
    FArrayBox yd;
    /** z coordinates of nodes */
    FArrayBox zd;

    /** normals of the x faces */
    FArrayBox xface_normals;
    /** normals of the y faces */
    FArrayBox yface_normals;
    /** tangents of the x faces */
    FArrayBox xface_tangents;
    /** tangents of the y faces */
    FArrayBox yface_tangents;
    /** normals of cell surfaces */
    FArrayBox surf_normals;
    /** edge lengths for both x and y faces */
    FArrayBox edge_lengths;

    /** area of each cell */
    FArrayBox area;
    /** curvature of each cell */
    FArrayBox curvature;
};

/**
 * @brief Allocate a new metric patch
 * 
 * @return fclaw2d_metric_patch_t* the new metric patch
 */
fclaw2d_metric_patch_t* fclaw2d_metric_patch_new();

/**
 * @brief Delete the metric patch
 * 
 * @param patchmp the metric patch, null on return
 */
void fclaw2d_metric_patch_delete(fclaw2d_metric_patch_t **patchmp);

#endif /* !FCLAW2D_METRIC_HPP */
