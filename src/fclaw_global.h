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

#ifndef FCLAW_GLOBAL_H
#define FCLAW_GLOBAL_H

#include <forestclaw.h>
#include <fclaw_timer.h>   /* Needed to create statically allocated array of timers */
#include <fclaw_pointer_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct fclaw_global fclaw_global_t;
typedef struct fclaw2d_map_context fclaw2d_map_context_t;
typedef struct fclaw3d_map_context fclaw3d_map_context_t;

struct fclaw_global
{
    int count_amr_advance;
    int count_ghost_exchange;
    int count_amr_regrid;
    int count_amr_new_domain;
    int count_single_step;
    int count_elliptic_grids;
    int count_multiproc_corner;
    int count_grids_per_proc;
    int count_grids_remote_boundary;
    int count_grids_local_boundary;
    fclaw_timer_t timers[FCLAW_TIMER_COUNT];

    /* Time at start of each subcycled time step */
    double curr_time;
    double curr_dt;

    sc_MPI_Comm mpicomm;
    int mpisize;              /**< Size of communicator. */
    int mpirank;              /**< Rank of this process in \b mpicomm. */

    /** Solver packages for internal use. */
    struct fclaw_package_container *pkg_container;

    struct fclaw_pointer_map *vtables;    /**< Vtables */
    struct fclaw_pointer_map *options;    /**< options */
    struct fclaw_pointer_map *attributes;    /**< attributes, things that are not vtables, or options */

    struct fclaw_domain *domain;

    void *user;
};

typedef struct fclaw_global_iterate fclaw_global_iterate_t;

struct fclaw_global_iterate
{
    fclaw_global_t* glob;
    void* user;
};

/** Allocate a new global structure. */
fclaw_global_t* fclaw_global_new (void);

fclaw_global_t* fclaw_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank);

void fclaw_global_destroy (fclaw_global_t * glob);

void fclaw_global_store_domain (fclaw_global_t* glob,
                                  struct fclaw_domain* domain);

void fclaw_global_iterate_level (fclaw_global_t * glob, int level,
                                   fclaw_patch_callback_t pcb, void *user);

void fclaw_global_iterate_patches (fclaw_global_t * glob,
                                     fclaw_patch_callback_t pcb, void *user);

void fclaw_global_iterate_families (fclaw_global_t * glob,
                                      fclaw_patch_callback_t pcb, void *user);

void fclaw_global_iterate_adapted (fclaw_global_t * glob,
                                     struct fclaw_domain* new_domain,
                                     fclaw_match_callback_t mcb, void *user);

void fclaw_global_iterate_level_mthread (fclaw_global_t * glob, int level,
                                           fclaw_patch_callback_t pcb, void *user);

void fclaw_global_iterate_partitioned (fclaw_global_t * glob,
                                         struct fclaw_domain * new_domain,
                                         fclaw_transfer_callback_t tcb,
                                         void *user);
/**
 * @brief Pack global structure into buffer
 * 
 * @param glob the global structure
 * @param buffer the buffer to write to
 * @return size_t number of bytes written
 */
size_t fclaw_global_pack(const fclaw_global_t * glob, char* buffer);

/**
 * @brief Get the number of bytes needed to pack the global structure
 * 
 * @param glob the structure
 * @return size_t the number of bytes needed to store structure
 */
size_t fclaw_global_packsize(const fclaw_global_t * glob);

/**
 * @brief Unpack global structure from buffer
 * 
 * @param buffer the buffer to read from
 * @param glob newly create global structure
 * @return size_t number of bytes read
 */
size_t fclaw_global_unpack(char* buffer, fclaw_global_t** glob);

/**
 * @brief Store an options structure in the glob
 * 
 * @param glob the global context
 * @param key the key to store the options under
 * @param options the options structure
 */
void fclaw_global_options_store (fclaw_global_t* glob, const char* key, void* options);

/**
 * @brief Get an options structure from the glob
 * 
 * @param glob the global context
 * @param key the key to retrieve the options from
 * @return void* the options
 */
void* fclaw_global_get_options (fclaw_global_t* glob, const char* key);

/**
 * @brief Store an attribute in the glob
 * 
 * @param glob the global context
 * @param key the key to store the attribute under
 * @param attrubute the attribute to store
 * @param destory callback to destroy the attribute. Optional, can be set to NULL
 */
void fclaw_global_attribute_store (fclaw_global_t* glob, 
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
void* fclaw_global_get_attribute (fclaw_global_t* glob, const char* key);

/**
 * @brief Store a glob variable in static memory
 *
 * @param glob the glob variable
 */
void fclaw_global_set_static (fclaw_global_t* glob);

/**
 * @brief Set the static glob variable to NULL
 */
void fclaw_global_clear_static (void);

/**
 * @brief Get the static glob variable
 *
 * @return fclaw_global_t* the glob variable
 */
fclaw_global_t* fclaw_global_get_static_global (void);

/**
 * @brief

 * @param glob
 */
void fclaw_set_global_context(fclaw_global_t *glob);

/**
 * @brief
 *
 * @param glob
 */
void fclaw_clear_global_context(fclaw_global_t *glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW_GLOBAL_H */
