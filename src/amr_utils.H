/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#ifndef AMR_UTILS_H
#define AMR_UTILS_H

#include "fclaw2d_map.H"
#include "fclaw_options.h"
#include "forestclaw2d.h"

#include "fclaw2d_typedefs.h"

#include "fclaw2d_capi.h"

/* Note: either we make this a C .h file, or we remove the extern "C". */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


// -----------------------------------------------------------------
// Check parms
// -----------------------------------------------------------------

void check_amr_parms (amr_options_t * gparms);

/* -----------------------------------------------------------------
   Initialize data
   ----------------------------------------------------------------- */
void init_domain_data (fclaw2d_domain_t * domain);
void init_block_data (fclaw2d_block_t * block);
void init_patch_data (fclaw2d_patch_t * patch);
void init_block_and_patch_data (fclaw2d_domain_t * domain);

void delete_domain_data (fclaw2d_domain_t * domain);
void delete_patch_data (fclaw2d_patch_t * patch);

/* -----------------------------------------------------------------
   Work with timers
   ----------------------------------------------------------------- */

/** Get current wall clock time in seconds.
 */
double fclaw2d_timer_wtime (void);
void fclaw2d_timer_init (fclaw2d_timer_t * timer);
void fclaw2d_timer_start (fclaw2d_timer_t * timer);
void fclaw2d_timer_stop (fclaw2d_timer_t * timer);

/* -----------------------------------------------------------------
   Return pointer to user data
   ---------------------------------------------------------------- */
fclaw2d_domain_data_t *get_domain_data (fclaw2d_domain_t * domain);
fclaw2d_block_data_t *get_block_data (fclaw2d_block_t * block);
fclaw2d_patch_data_t *get_patch_data (fclaw2d_patch_t * patch);
fclaw2d_map_context_t *get_map_context(fclaw2d_domain_t * domain);

int fclaw2d_domain_get_num_patches(fclaw2d_domain_t* domain);


/* -----------------------------------------------------------------
   Set user data with user defined variables, etc.
   ----------------------------------------------------------------- */
void set_block_data (fclaw2d_block_t * block, const int mthbc[]);

void copy_domain_data (fclaw2d_domain_t * old_domain,
                       fclaw2d_domain_t * new_domain);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
