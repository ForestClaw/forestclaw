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

typedef struct fclaw_patch_callback_wrap_user
{
    fclaw_patch_callback_t pcb;
    void* user;
} fclaw_patch_callback_wrap_user_t;

typedef struct fclaw_transfer_callback_wrap_user
{
    fclaw_transfer_callback_t tcb;
    void *user;
} fclaw_transfer_callback_wrap_user_t;

typedef struct fclaw_match_callback_wrap_user
{
    fclaw_match_callback_t mcb;
    void *user;
} fclaw_match_callback_wrap_user_t;

typedef struct fclaw_integrate_ray_wrap_user
{
    fclaw_integrate_ray_t intersect;
    void *user;
} fclaw_integrate_ray_wrap_user_t;

typedef struct fclaw_interpolate_pointer_wrap_user
{
    fclaw_interpolate_point_t interpolate;
    void *user;
} fclaw_interpolate_point_user_wrap_t;

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
