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
 
    struct fclaw_package_container *pkg_container;    /**< Solver packages for internal use. */
 
    struct fclaw2d_map_context* cont;
    struct fclaw2d_domain *domain;

    struct fclaw2d_diagnostics_accumulator *acc;

    void *user;
};

struct fclaw2d_global_iterate
{
    fclaw2d_global_t* glob;
    void* user;

};

/* Use forward references here, since this file gets included everywhere */
struct fclaw2d_domain;
struct fclaw2d_map_context;
struct fclaw_package_container;
struct fclaw2d_diagnostics_accumulator;

/** Allocate a new global structure.
 * \param [in] gparms           If not NULL, we borrow this gparms pointer.
 *                              If NULL, we allocate gparms ourselves.
 */

fclaw2d_global_t* fclaw2d_global_new (void);

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

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_GLOBAL_H */
