/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW3D_GLOBAL_H
#define FCLAW3D_GLOBAL_H

#include <fclaw_global.h>
#include <forestclaw3d.h>  /* Needed to declare callbacks (below) */
#include <fclaw3d_map.h>   /* Needed to store the map context */

#include <fclaw_timer.h>   /* Needed to create statically allocated array of timers */
#include <fclaw_pointer_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* these are dimension-specific functions */

void fclaw3d_iterate_patch_cb
    (fclaw_domain_t * domain, fclaw_patch_t * patch,
     int blockno, int patchno, void *user);

void fclaw3d_iterate_family_cb
    (fclaw_domain_t * domain, fclaw_patch_t * patch,
     int blockno, int patchno, void *user);

/* much of the following will move into fclaw_global.h */

void fclaw3d_global_store_map (fclaw_global_t* glob,
                               fclaw_map_context_t * map);

fclaw_map_context_t* fclaw3d_global_get_map(fclaw_global_t* glob);

/**
 * @brief Set the static glob variable to NULL
 */
void fclaw3d_global_unset_global (void);

/**
 * @brief Get the static glob variable
 *
 * @return fclaw_global_t* the glob variable
 */
fclaw_global_t* fclaw3d_global_get_global (void);

/**
 * @brief
 *
 * @param glob
 */
void fclaw3d_set_global_context(fclaw_global_t *glob);

/**
 * @brief
 *
 * @param glob
 */
void fclaw3d_clear_global_context(fclaw_global_t *glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW3D_GLOBAL_H */
