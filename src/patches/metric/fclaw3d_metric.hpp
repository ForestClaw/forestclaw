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

#ifndef FCLAW3D_METRIC_HPP
#define FCLAW3D_METRIC_HPP

#include <fclaw_farraybox.hpp>     /* Needed for FArray box used to store metric data */

/**
 * @brief Struct for patch metric data
 */
struct fclaw3d_metric_patch_t
{
    /** The number of cells in the x direction */
    int mx;           
    /** The number of cells in the x direction */
    int my;           
    /** The number of cells in the z direction */
    int mz;           

    /** The number of ghost cells */
    int mbc;
    /** The block number */
    int blockno;

    /** The spacing in the x direction */
    double dx;
    /** The spacing in the y direction */
    double dy;
    /** The spacing in the z direction */
    double dz;

    /** x value of corner (0,0,0) in 3d patch */
    double xlower;
    /** y value of corner (0,0,0) in 3d patch */
    double ylower;
    /** z value of corner (0,0,0) in 3d patch */
    double zlower;

    /** x value of corner (1,1,1) in 3d patch */
    double xupper;
    /** y value of corner (1,1,1) in 3d patch */
    double yupper;
    /** z value of corner (1,1,1) in 3d patch */
    double zupper; 

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

    /** Rotation matrix at face with constant x (left) */
    FArrayBox xrot;

    /** Rotation matrix at face with constant y (front) */
    FArrayBox yrot;

    /** Rotation matrix at face with constant z (bottom) */
    FArrayBox zrot;

    /** volume of each cell */
    FArrayBox volume;

    /** Area of each cell of three faces */
    FArrayBox face_area;

};

/**
 * @brief Allocate a new metric patch
 * 
 * @return fclaw2d_metric_patch_t* the new metric patch
 */
fclaw3d_metric_patch_t* fclaw3d_metric_patch_new();

/**
 * @brief Delete the metric patch
 * 
 * @param patchmp the metric patch, null on return
 */
void fclaw3d_metric_patch_delete(fclaw3d_metric_patch_t **patchmp);

#endif /* !FCLAW3D_METRIC_HPP */
