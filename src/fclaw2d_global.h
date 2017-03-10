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

#include <forestclaw2d.h>
#include <fclaw_options.h>
#include <fclaw_package.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

#define FCLAW2D_SPACEDIM 2
extern const int SpaceDim;

/* Number of faces to a patch. Changed from CUBEFACES to NUMFACES to
   avoid any confusion in the 2d case. */
#define FCLAW2D_NUMFACES (2 * FCLAW2D_SPACEDIM)
extern const int NumFaces;

#define FCLAW2D_P4EST_REFINE_FACTOR 2
extern const int p4est_refineFactor;

#define FCLAW2D_NUM_CORNERS 4
extern const int NumCorners;

#define FCLAW2D_NUM_SIBLINGS 4
extern const int NumSiblings;

typedef struct fclaw2d_global
{

#if 0
    int count_set_clawpatch, count_delete_clawpatch;
    int count_amr_advance, count_ghost_exchange, count_amr_regrid;
    int count_amr_new_domain;
    int count_single_step;
    int count_multiproc_corner;
    int count_grids_per_proc;
    int count_grids_remote_boundary;
    int count_grids_local_boundary;
    int is_latest_domain;
    fclaw2d_timer_t timers[FCLAW2D_TIMER_COUNT];

    /* Time at start of each subcycled time step */
    double curr_time;
#endif
	sc_MPI_Comm mpicomm;
    int mpisize;              /**< Size of communicator. */
    int mpirank;              /**< Rank of this process in \b mpicomm. */

    int gparms_owned;                   /**< Did we allocate \a gparms? */
    fclaw_options_t *gparms;            /**< Option values for forestclaw. */
 
    fclaw_package_container_t *pkgs;    /**< Solver packages for internal use. */
 
    fclaw2d_domain_t *domain;
}
fclaw2d_global_t;

typedef struct fclaw2d_global_iterate
{
    fclaw2d_global_t* glob;
    void* user;
}fclaw2d_global_iterate_t;

/** Allocate a new global structure.
 * \param [in] gparms           If not NULL, we borrow this gparms pointer.
 *                              If NULL, we allocate gparms ourselves.
 */
fclaw2d_global_t *fclaw2d_global_new (fclaw_options_t * gparms, fclaw2d_domain_t* domain);

/** Free a global structures and all members. */
void fclaw2d_global_destroy (fclaw2d_global_t * glob);

/** Access the package container from the global type. */
fclaw_package_container_t *fclaw2d_global_get_container (fclaw2d_global_t *
                                                         glob);

fclaw_app_t* fclaw2d_global_get_app(fclaw2d_global_t* glob);

void fclaw2d_global_iterate_level (fclaw2d_global_t * glob, int level,
                                   fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_patches (fclaw2d_global_t * glob,
                                     fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_families (fclaw2d_global_t * glob,
                                      fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_global_iterate_adapted (fclaw2d_global_t * glob, fclaw2d_domain_t* new_domain,
                                     fclaw2d_match_callback_t mcb, void *user);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_GLOBAL_H */
