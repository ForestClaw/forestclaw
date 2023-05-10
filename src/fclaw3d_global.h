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

#ifndef FCLAW3D_GLOBAL_H
#define FCLAW3D_GLOBAL_H

#include <forestclaw3d.h>  /* Needed to declare callbacks (below) */
#include <fclaw3d_map.h>   /* Needed to store the map context */

#include <fclaw_timer.h>   /* Needed to create statically allocated array of timers */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* these are dimension-specific functions */

void fclaw3d_iterate_patch_cb
    (fclaw3d_domain_t * domain, fclaw3d_patch_t * patch,
     int blockno, int patchno, void *user);

void fclaw3d_iterate_family_cb
    (fclaw3d_domain_t * domain, fclaw3d_patch_t * patch,
     int blockno, int patchno, void *user);

/* much of the following will move into fclaw_global.h */

typedef struct fclaw3d_global fclaw3d_global_t;
typedef struct fclaw3d_global_iterate fclaw3d_global_iterate_t;

struct fclaw3d_global
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

    struct fclaw3d_map_context* cont;
    struct fclaw3d_domain *domain;

#if 0
    /* CB: is this a good place for the accumulator?
           Would it be possible to add and retrieve it as an anonymous
           object that does not need to be known to this file? */

    struct fclaw3d_diagnostics_accumulator *acc;
#endif

    void *user;
};

struct fclaw3d_global_iterate
{
    fclaw3d_global_t* glob;
    void* user;
};

/** Allocate a new global structure. */
fclaw3d_global_t* fclaw3d_global_new (void);

fclaw3d_global_t* fclaw3d_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank);

void fclaw3d_global_destroy (fclaw3d_global_t * glob);

void fclaw3d_global_store_domain (fclaw3d_global_t* glob,
                                  struct fclaw3d_domain* domain);

void fclaw3d_global_store_map (fclaw3d_global_t* glob,
                               fclaw3d_map_context_t * map);

void fclaw3d_global_iterate_level (fclaw3d_global_t * glob, int level,
                                   fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_patches (fclaw3d_global_t * glob,
                                     fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_families (fclaw3d_global_t * glob,
                                      fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_adapted (fclaw3d_global_t * glob,
                                     struct fclaw3d_domain* new_domain,
                                     fclaw3d_match_callback_t mcb, void *user);

void fclaw3d_global_iterate_level_mthread (fclaw3d_global_t * glob, int level,
                                           fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_partitioned (fclaw3d_global_t * glob,
                                         struct fclaw3d_domain * new_domain,
                                         fclaw3d_transfer_callback_t tcb,
                                         void *user);

/**
 * @brief Store a glob variable in static memory
 *
 * @param glob the glob variable
 */
void fclaw3d_global_set_global (fclaw3d_global_t* glob);

/**
 * @brief Set the static glob variable to NULL
 */
void fclaw3d_global_unset_global (void);

/**
 * @brief Get the static glob variable
 *
 * @return fclaw2d_global_t* the glob variable
 */
fclaw3d_global_t* fclaw3d_global_get_global (void);

/**
 * @brief
 *
 * @param glob
 */
void fclaw3d_set_global_context(fclaw3d_global_t *glob);

/**
 * @brief
 *
 * @param glob
 */
void fclaw3d_clear_global_context(fclaw3d_global_t *glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW3D_GLOBAL_H */
