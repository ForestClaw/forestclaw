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

#ifndef FCLAW2D_GHOST_FILL_H
#define FCLAW2D_GHOST_FILL_H

#include <fclaw_timer.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_domain;
struct fclaw2d_patch;

/* -------------------------------------------
   Routines needed to fill in ghost cells
   ------------------------------------------- */

typedef enum fclaw2d_ghost_fill_parallel_mode
{
	FCLAW2D_BOUNDARY_INTERIOR_ONLY = 0,  /* Don't read parallel patches */
	FCLAW2D_BOUNDARY_GHOST_ONLY,         /* read parallel patches */
	FCLAW2D_BOUNDARY_ALL,                 /* read parallel patches */
	FCLAW2D_BOUNDARY_LOCAL_ALL,
} fclaw2d_ghost_fill_parallel_mode_t;


typedef enum fclaw2d_exchange_type
{
	FCLAW2D_COPY = 1,
	FCLAW2D_AVERAGE,
	FCLAW2D_INTERPOLATE,
	FCLAW2D_TIME_SYNC_F2C,
	FCLAW2D_TIME_SYNC_SAMESIZE,
} fclaw2d_exchange_type_t;

typedef enum fclaw2d_grid_type
{
	FCLAW2D_IS_COARSE = 1,
	FCLAW2D_IS_FINE,
} fclaw2d_grid_type_t;

typedef struct fclaw2d_exchange_info
{
	int time_interp;
	int level;
	int read_parallel_patches;   /* before we have done a parallel exchange */
	fclaw2d_exchange_type_t exchange_type;
	fclaw2d_grid_type_t grid_type;
	int has_fine_grid_neighbor;

} fclaw2d_exchange_info_t;

void cb_corner_fill(struct fclaw2d_domain *domain,
					struct fclaw2d_patch *this_patch,
					int this_block_idx,
					int this_patch_idx,
					void *user);

void cb_face_fill(struct fclaw2d_domain *domain,
				  struct fclaw2d_patch *this_patch,
				  int this_block_idx,
				  int this_patch_idx,
				  void *user);

void fclaw2d_ghost_update(struct fclaw2d_global* glob,
						  int fine_level,
						  int coarse_level,
						  double sync_time,
						  int time_interp,
						  fclaw2d_timer_names_t running);

void fclaw2d_ghost_update_async(struct fclaw2d_global* glob,
								int fine_level,
								int coarse_level,
								double sync_time,
								int time_interp,
								fclaw2d_timer_names_t running);

void fclaw2d_ghost_update_nonasync(struct fclaw2d_global* glob,
								   int fine_level,
								   int coarse_level,
								   double sync_time,
								   int time_interp,
								   fclaw2d_timer_names_t running);

void fclaw2d_face_neighbor_ghost(struct fclaw2d_global* glob,
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
