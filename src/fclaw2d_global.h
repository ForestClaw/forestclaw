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

#ifndef FCLAW2D_GLOBAL_H
#define FCLAW2D_GLOBAL_H

#include <forestclaw2d.h>  /* Needed to declare callbacks (below) */
#include <fclaw_global.h>
#include <fclaw2d_map.h>   /* Needed to store the map context */

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

void fclaw2d_iterate_patch_cb
  (fclaw_domain_t *domain, fclaw_patch_t *patch,
   int blockno, int patchno, void *user);

void fclaw2d_iterate_family_cb
  (fclaw_domain_t *domain, fclaw_patch_t *patch,
   int blockno, int patchno, void *user);

/* much of the following will move into fclaw_global.h */



/* Use forward references here, since this file gets included everywhere */
/* CB: is there a way not to need the forward references?
       Depending on fclaw2d_domain, _map, package seems entirely acceptable.
       For those, including the respective headers might not be so bad.
       About the diagnostics accumulator see remark above. */
struct fclaw_domain;
struct fclaw2d_map_context;
struct fclaw_package_container;
struct fclaw2d_diagnostics_accumulator;

void fclaw2d_global_store_map (fclaw_global_t* glob,
                               struct fclaw2d_map_context * map);

fclaw2d_map_context_t* fclaw2d_global_get_map(fclaw_global_t* glob);

/**
 * @brief

 * @param glob
 */
void fclaw2d_set_global_context(fclaw_global_t *glob);

/**
 * @brief
 *
 * @param glob
 */
void fclaw2d_clear_global_context(fclaw_global_t *glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_GLOBAL_H */
