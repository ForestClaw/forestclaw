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

#include <fclaw2d_time_sync.h>

#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_options.h>
#include <fclaw2d_exchange.h>

void fclaw2d_time_sync_reset(fclaw2d_global_t *glob, 
							 int minlevel,int maxlevel, int init)
{
	/* This is used for updating conservation arrays, for example */
	fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();

	// FCLAW_ASSERT(minlevel == maxlevel);
	// fclaw_global_essentialf("reset level = %d\n",minlevel);
	if (fclaw_vt->time_sync_reset != NULL)
	{
		fclaw_vt->time_sync_reset(glob,minlevel,maxlevel,init);
	}
}


static
void copy_at_blockbdry(fclaw2d_global_t *glob,
					   int level,
					   int read_parallel_patches,
					   fclaw2d_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw2d_exchange_info_t e_info;
	e_info.exchange_type = FCLAW2D_TIME_SYNC_COPY;
	e_info.grid_type = FCLAW2D_IS_COARSE;
	e_info.time_interp = 0;
	e_info.read_parallel_patches = read_parallel_patches;

	fclaw2d_global_iterate_level(glob, level, cb_face_fill,
						 &e_info);    
}


static
void fine2coarse(fclaw2d_global_t *glob,
				int level,
				int read_parallel_patches,
				fclaw2d_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw2d_exchange_info_t e_info;
	e_info.exchange_type = FCLAW2D_TIME_SYNC_FINE_TO_COARSE;
	e_info.grid_type = FCLAW2D_IS_COARSE;
	e_info.time_interp = 0;
	e_info.read_parallel_patches = read_parallel_patches;

	fclaw2d_global_iterate_level(glob, level, cb_face_fill,
	                                &e_info);
}

static
void correct_coarse_cells(fclaw2d_global_t *glob, 
                          int minlevel, 
                          int read_parallel_patches,
                          fclaw2d_ghost_fill_parallel_mode_t ghost_mode)

{
	fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

	int level;

	int mincoarse = minlevel;
	int maxcoarse = fclaw_opt->maxlevel - 1;
	
	/* We have to at all coarse/fine boundaries, and at block_boundaries
	where grids are the same size */
	for(level = maxcoarse; level >= mincoarse; level--)
	{
		/* Level is the coarse level to be corrected */
		fine2coarse(glob,level,
					read_parallel_patches,ghost_mode);

		/* Reset fluxes on fine grid, since they have all been used to 
		   correct coarse grid.  Coarse grid data is reset in
		   syncing step. */
		fclaw2d_time_sync_reset(glob, level+1, level+1, 1);
	}

#if 0
	/* This step accounts for any metric discontinuities at block boundaries */
	for(level = maxcoarse; level >= mincoarse; level--)
	{
		copy_at_blockbdry(glob,level,
						  read_parallel_patches,ghost_mode);
	}
#endif	

}


void fclaw2d_time_sync(fclaw2d_global_t *glob, int minlevel, int maxlevel)
{

	int time_interp = 0;

	fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_TIMESYNC]);


#if 0
	/* --------------------------------------------------------------
		Send and receive patches
	-----------------------------------------------------------------*/
	fclaw2d_exchange_ghost_patches_begin(glob,minlevel,maxlevel,time_interp,
										 FCLAW2D_TIMER_TIMESYNC);

	fclaw2d_exchange_ghost_patches_end(glob,minlevel,maxlevel,time_interp,
									   FCLAW2D_TIMER_TIMESYNC);

	/* Three-way corner exchanges */
	fclaw2d_face_neighbor_ghost(glob,minlevel,maxlevel,time_interp);
#endif

	/* -------------------------------------------------------------
	Add corrections from  fine grids to coarse grid.  This is is done 
	when two or more levels are time synchronized.
	   -- All parallel patches are valid 
	   -- Iterate over boundary patches only, since correction occurs
	   only at coarse/fine boundaries.
	------------------------------------------------------------- */

	/* Parallel patches are valid */
	int read_parallel_patches = 1;

	/* Indicates that we should read only interior, only boundary, or all patches
	on a local processor */
	fclaw2d_ghost_fill_parallel_mode_t parallel_mode =
		   FCLAW2D_BOUNDARY_ALL;    

	correct_coarse_cells(glob,minlevel,read_parallel_patches,parallel_mode);


	fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_TIMESYNC]);

}


