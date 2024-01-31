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

#ifndef FCLAW_GHOST_FILL_H
#define FCLAW_GHOST_FILL_H

#include <fclaw_timer.h>    /* Needed for definition of timer names */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw_global;
struct fclaw_domain;
struct fclaw_patch;

/** 
 *  @file
 *  Routines needed to fill in ghost cells
 */

typedef enum fclaw_ghost_fill_parallel_mode
{
	FCLAW_BOUNDARY_INTERIOR_ONLY = 0,  /* Don't read parallel patches */
	FCLAW_BOUNDARY_GHOST_ONLY,         /* read parallel patches */
	FCLAW_BOUNDARY_ALL,                 /* read parallel patches */
	FCLAW_BOUNDARY_LOCAL_ALL,
} fclaw_ghost_fill_parallel_mode_t;


typedef enum fclaw_exchange_type
{
	FCLAW_COPY = 1,
	FCLAW_AVERAGE,
	FCLAW_INTERPOLATE,
	FCLAW_TIME_SYNC_F2C,
	FCLAW_TIME_SYNC_SAMESIZE,
} fclaw_exchange_type_t;

/** enum for which type of grid is being filled */
typedef enum fclaw_grid_type
{
	/** fill ghost cells for coarser neighbor grid */
	FCLAW_IS_COARSE = 1,
	/** fill ghost cells for finer neighbor grid */
	FCLAW_IS_FINE,
} fclaw_grid_type_t;

typedef struct fclaw_exchange_info
{
	int time_interp;
	int level;
	int read_parallel_patches;   /* before we have done a parallel exchange */
	fclaw_exchange_type_t exchange_type;
	fclaw_grid_type_t grid_type;
	int has_fine_grid_neighbor;

} fclaw_exchange_info_t;

void cb_corner_fill(struct fclaw_domain *domain,
					struct fclaw_patch *this_patch,
					int this_block_idx,
					int this_patch_idx,
					void *user);

void cb_face_fill(struct fclaw_domain *domain,
				  struct fclaw_patch *this_patch,
				  int this_block_idx,
				  int this_patch_idx,
				  void *user);

void fclaw_ghost_update(struct fclaw_global* glob,
						  int fine_level,
						  int coarse_level,
						  double sync_time,
						  int time_interp,
						  fclaw_timer_names_t running);

void fclaw_ghost_update_async(struct fclaw_global* glob,
								int fine_level,
								int coarse_level,
								double sync_time,
								int time_interp,
								fclaw_timer_names_t running);

/**
 * <summary>Complete exchange of all ghost patches at all levels.</summary>
 * <remarks>All parallel ghost patches are also exchanged at all
 * levels.</remarks>
 * <list>
 *    <item>Every level exchanges ghost cells with other patches
 *       at that level</item>
 *    <item>Every finer level exchanges with a coarser level</item>
 *    <item> All levels will be updated in next update step, regardless of
 *       whether we are in the subcycled or non-subcycled case.</item>
 *       </list>
 *   The reason for two separate ghost cell exchange routines is that
 *   the logic here is considerably simpler than for the partial
 *   update used in intermediate steps in the subcycled case.
 **/
void fclaw_ghost_update_nonasync(struct fclaw_global* glob,
								   int fine_level,
								   int coarse_level,
								   double sync_time,
								   int time_interp,
								   fclaw_timer_names_t running);

/**
 * <summary>Complete exchange of all ghost patches at all levels.</summary>
 * <remarks>All parallel ghost patches are also exchanged at all
 * levels.</remarks>
 * <list>
 *    <item>Every level exchanges ghost cells with other patches
 *       at that level</item>
 *    <item>Every finer level exchanges with a coarser level</item>
 *    <item> All levels will be updated in next update step, regardless of
 *       whether we are in the subcycled or non-subcycled case.</item>
 *       </list>
 *   The reason for two separate ghost cell exchange routines is that
 *   the logic here is considerably simpler than for the partial
 *   update used in intermediate steps in the subcycled case.
 **/
void fclaw_face_neighbor_ghost(struct fclaw_global* glob,
								 int minlevel,
								 int maxlevel,
								 int time_interp);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
