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

#ifndef FCLAW2D_GLOBAL_H
#define FCLAW2D_GLOBAL_H

#include <forestclaw2d.h>  /* Needed to declare callbacks (below) */
#include <fclaw2d_map.h>   /* Needed to store the map context */

#include <fclaw_timer.h>   /* Needed to create statically allocated array of timers */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* these are dimension-specific functions */

void fclaw2d_iterate_patch_cb
  (fclaw2d_domain_t *domain, fclaw2d_patch_t *patch,
   int blockno, int patchno, void *user);

void fclaw2d_iterate_family_cb
  (fclaw2d_domain_t *domain, fclaw2d_patch_t *patch,
   int blockno, int patchno, void *user);

/* much of the following will move into fclaw_global.h */

typedef struct fclaw2d_global fclaw2d_global_t;
typedef struct fclaw2d_global_iterate fclaw2d_global_iterate_t;

struct fclaw2d_global
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
    fclaw2d_timer_t timers[FCLAW2D_TIMER_COUNT];

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

    struct fclaw2d_map_context* cont;
    struct fclaw2d_domain *domain;

    /* CB: is this a good place for the accumulator?
           Would it be possible to add and retrieve it as an anonymous
           object that does not need to be known to this file? */
    struct fclaw2d_diagnostics_accumulator *acc;

    /* CB: this is application specific.
           Would it not be cleaner to add the gauges in a way to global
           that this file does not need to know about gauges at all? */
    struct fclaw_gauge_info* gauge_info;

    void *user;
};

struct fclaw2d_global_iterate
{
    fclaw2d_global_t* glob;
    void* user;
};

/* Use forward references here, since this file gets included everywhere */
/* CB: is there a way not to need the forward references?
       Depending on fclaw2d_domain, _map, package seems entirely acceptable.
       For those, including the respective headers might not be so bad.
       About the diagnostics accumulator see remark above. */
struct fclaw2d_domain;
struct fclaw2d_map_context;
struct fclaw_package_container;
struct fclaw2d_diagnostics_accumulator;

/** Allocate a new global structure. */
fclaw2d_global_t* fclaw2d_global_new (void);

fclaw2d_global_t* fclaw2d_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank);

void fclaw2d_global_destroy (fclaw2d_global_t * glob);

void fclaw2d_global_store_domain (fclaw2d_global_t* glob,
                                  struct fclaw2d_domain* domain);

void fclaw2d_global_store_map (fclaw2d_global_t* glob,
                               struct fclaw2d_map_context * map);

void fclaw2d_global_iterate_level (fclaw2d_global_t * glob, int level,
                                   fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_patches (fclaw2d_global_t * glob,
                                     fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_families (fclaw2d_global_t * glob,
                                      fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_adapted (fclaw2d_global_t * glob,
                                     struct fclaw2d_domain* new_domain,
                                     fclaw2d_match_callback_t mcb, void *user);

void fclaw2d_global_iterate_level_mthread (fclaw2d_global_t * glob, int level,
                                           fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_partitioned (fclaw2d_global_t * glob,
                                         struct fclaw2d_domain * new_domain,
                                         fclaw2d_transfer_callback_t tcb,
                                         void *user);

/**
 * @brief Store a glob variable in static memory
 *
 * @param glob the glob variable
 */
void fclaw2d_global_set_global (fclaw2d_global_t* glob);

/**
 * @brief Set the static glob variable to NULL
 */
void fclaw2d_global_unset_global (void);

/**
 * @brief Get the static glob variable
 *
 * @return fclaw2d_global_t* the glob variable
 */
fclaw2d_global_t* fclaw2d_global_get_global (void);

/**
 * @brief
 *
 * @param glob
 */
void fclaw2d_set_global_context(fclaw2d_global_t *glob);

/**
 * @brief
 *
 * @param glob
 */
void fclaw2d_clear_global_context(fclaw2d_global_t *glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_GLOBAL_H */
