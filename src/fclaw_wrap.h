/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
 * @brief Types needed for underlying wrapped types
 *        Most users will not need to include the file, unless doing something more advanced.
 */


#ifndef FCLAW_WRAP_H
#define FCLAW_WRAP_H

#include <forestclaw.h>
#include <fclaw_convenience.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/**
 * @brief Callback wrapper struct for patch callbacks.
 *
 * This allows callbacks with dimension independent types to be 
 * called from dimension dependent code.
 *
 * This type should be passed in as the user pointer alongside
 * @ref fclaw2d_patch_wrap_callback or @ref fclaw3d_patch_wrap_callback
 * to the function that takes a dimensioned callback.
 *
 * @ref fclaw2d_patch_wrap_callback or @ref fclaw3d_patch_wrap_callback
 * will then call the dimension independent callback specified in the struct,
 * passing the user pointer specified in the struct.
 * 
 */
typedef struct fclaw_patch_wrap_user
{
    /** Dimension independent patch callback to call */
    fclaw_patch_callback_t pcb;
    /** User pointer to pass to dimension independent callback */
    void* user;
} fclaw_patch_wrap_user_t;

/**
 * @brief Callback wrapper struct for transfer callbacks.
 *
 * This allows callbacks with dimension independent types to be 
 * called from dimension dependent code.
 *
 * This type should be passed in as the user pointer alongside
 * @ref fclaw2d_transfer_wrap_callback or @ref fclaw3d_transfer_wrap_callback
 * to the function that takes a dimensioned callback.
 *
 * @ref fclaw2d_transfer_wrap_callback or @ref fclaw3d_transfer_wrap_callback
 * will then call the dimension independent callback specified in the struct,
 * passing the user pointer specified in the struct.
 * 
 */
typedef struct fclaw_transfer_wrap_user
{
    /** Dimension independent transfer callback to call */
    fclaw_transfer_callback_t tcb;
    /** User pointer to pass to dimension independent callback */
    void *user;
} fclaw_transfer_wrap_user_t;

/**
 * @brief Callback wrapper struct for match callbacks.
 *
 * This allows callbacks with dimension independent types to be 
 * called from dimension dependent code.
 *
 * This type should be passed in as the user pointer alongside
 * @ref fclaw2d_match_wrap_callback or @ref fclaw3d_match_wrap_callback
 * to the function that takes a dimensioned callback.
 *
 * @ref fclaw2d_match_wrap_callback or @ref fclaw3d_match_wrap_callback
 * will then call the dimension independent callback specified in the struct,
 * passing the user pointer specified in the struct.
 * 
 */
typedef struct fclaw_match_wrap_user
{
    /** Dimension independent match callback to call */
    fclaw_match_callback_t mcb;
    /** User pointer to pass to dimension independent callback */
    void *user;
} fclaw_match_wrap_user_t;

/**
 * @brief Callback wrapper struct for match callbacks.
 *
 * This allows callbacks with dimension independent types to be 
 * called from dimension dependent code.
 *
 * This type should be passed in as the user pointer alongside
 * @ref fclaw2d_integrate_ray_wrap_callback or @ref fclaw3d_integrate_ray_wrap_callback
 * to the function that takes a dimensioned callback.
 *
 * @ref fclaw2d_integrate_ray_wrap_callback or @ref fclaw3d_integrate_ray_wrap_callback
 * will then call the dimension independent callback specified in the struct,
 * passing the user pointer specified in the struct.
 * 
 */
typedef struct fclaw_integrate_ray_wrap_user
{
    /** Dimension independent integrate ray callback to call */
    fclaw_integrate_ray_t intersect;
    /** User pointer to pass to dimension independent callback */
    void *user;
} fclaw_integrate_ray_wrap_user_t;

/**
 * @brief Callback wrapper struct for match callbacks.
 *
 * This allows callbacks with dimension independent types to be 
 * called from dimension dependent code.
 *
 * This type should be passed in as the user pointer alongside
 * @ref fclaw2d_iterpolate_point_wrap_callback or @ref fclaw3d_iterpolate_point_wrap_callback
 * to the function that takes a dimensioned callback.
 *
 * @ref fclaw2d_iterpolate_point_wrap_callback or @ref fclaw3d_iterpolate_point_wrap_callback
 * will then call the dimension independent callback specified in the struct,
 * passing the user pointer specified in the struct.
 * 
 */
typedef struct fclaw_interpolate_point_wrap_user
{
    /** Dimension independent interpolate callback to call */
    fclaw_interpolate_point_t interpolate;
    /** User pointer to pass to dimension independent callback */
    void *user;
} fclaw_interpolate_point_user_wrap_t;

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
