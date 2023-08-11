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

/** Allocate a new global structure. */
fclaw_global_t* fclaw2d_global_new (void);

fclaw_global_t* fclaw2d_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank);

void fclaw2d_global_destroy (fclaw_global_t * glob);

void fclaw2d_global_store_domain (fclaw_global_t* glob,
                                  struct fclaw_domain* domain);

void fclaw2d_global_store_map (fclaw_global_t* glob,
                               struct fclaw2d_map_context * map);

fclaw2d_map_context_t* fclaw2d_global_get_map(fclaw_global_t* glob);

void fclaw2d_global_iterate_level (fclaw_global_t * glob, int level,
                                   fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_patches (fclaw_global_t * glob,
                                     fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_families (fclaw_global_t * glob,
                                      fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_adapted (fclaw_global_t * glob,
                                     struct fclaw_domain* new_domain,
                                     fclaw2d_match_callback_t mcb, void *user);

void fclaw2d_global_iterate_level_mthread (fclaw_global_t * glob, int level,
                                           fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_partitioned (fclaw_global_t * glob,
                                         struct fclaw_domain * new_domain,
                                         fclaw2d_transfer_callback_t tcb,
                                         void *user);

/**
 * @brief Store an options structure in the glob
 * 
 * @param glob the global context
 * @param key the key to store the options under
 * @param options the options structure
 */
void fclaw2d_global_options_store (fclaw_global_t* glob, const char* key, void* options);

/**
 * @brief Get an options structure from the glob
 * 
 * @param glob the global context
 * @param key the key to retrieve the options from
 * @return void* the options
 */
void* fclaw2d_global_get_options (fclaw_global_t* glob, const char* key);

/**
 * @brief Store an attribute in the glob
 * 
 * @param glob the global context
 * @param key the key to store the attribute under
 * @param attrubute the attribute to store
 * @param destory callback to destroy the attribute. Optional, can be set to NULL
 */
void fclaw2d_global_attribute_store (fclaw_global_t* glob, 
                                     const char* key, 
                                     void* attribute,
                                     fclaw_pointer_map_value_destroy_t destroy);

/**
 * @brief Get an attribute structure from the glob
 * 
 * @param glob the global context
 * @param key the key to retrieve the attribute from
 * @return void* the options
 */
void* fclaw2d_global_get_attribute (fclaw_global_t* glob, const char* key);


/**
 * @brief Store a glob variable in static memory
 *
 * @param glob the glob variable
 */
void fclaw2d_global_set_global (fclaw_global_t* glob);

/**
 * @brief Set the static glob variable to NULL
 */
void fclaw2d_global_unset_global (void);

/**
 * @brief Get the static glob variable
 *
 * @return fclaw_global_t* the glob variable
 */
fclaw_global_t* fclaw2d_global_get_global (void);

/**
 * @brief
 *
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
