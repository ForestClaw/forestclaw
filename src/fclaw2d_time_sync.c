/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_options.h>
#include <fclaw2d_exchange.h>

typedef struct fclaw2d_time_sync_info
{
    /** The type of reset to perform */
    fclaw2d_time_sync_type_t reset_mode;
    /** The type of reset to perform */
    int coarse_level;
    /** The type of reset to perform */
    int minlevel;
} fclaw2d_time_sync_info_t;




static 
void cb_time_sync_reset(fclaw2d_domain_t *domain,
                        fclaw_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        void *user)
{
	fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;
	fclaw2d_global_t *glob = (fclaw2d_global_t*) g->glob;

	fclaw2d_time_sync_info_t *tstype = (fclaw2d_time_sync_info_t*) g->user;

	int reset_mode = tstype->reset_mode;
	int coarse_level = tstype->coarse_level;

	fclaw2d_patch_time_sync_reset(glob,this_patch,coarse_level,reset_mode);

}

static
void time_sync_reset (fclaw2d_global_t* glob,
                      int coarse_level, int reset_mode)
{
	fclaw2d_time_sync_info_t ts_info;

	ts_info.reset_mode = (fclaw2d_time_sync_type_t) reset_mode;
	ts_info.coarse_level = coarse_level;

	fclaw2d_global_iterate_level(glob, coarse_level, 
	                             cb_time_sync_reset, &ts_info);

	if (reset_mode == FCLAW2D_TIME_SYNC_RESET_F2C)
	{
		fclaw2d_global_iterate_level(glob, coarse_level+1, 
		                             cb_time_sync_reset, &ts_info);
	}

}

static
void copy_at_blockbdry(fclaw2d_global_t *glob,
					   int level,
					   int read_parallel_patches,
					   fclaw2d_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw2d_exchange_info_t e_info;
	e_info.exchange_type = FCLAW2D_TIME_SYNC_SAMESIZE;
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
	e_info.exchange_type = FCLAW2D_TIME_SYNC_F2C;
	e_info.grid_type = FCLAW2D_IS_COARSE;
	e_info.time_interp = 0;
	e_info.read_parallel_patches = read_parallel_patches;

	fclaw2d_global_iterate_level(glob, level, cb_face_fill, &e_info);
}

static
void correct_coarse_cells(fclaw2d_global_t *glob, 
                          int minlevel, 
                          int read_parallel_patches,
                          fclaw2d_ghost_fill_parallel_mode_t ghost_mode)

{
	fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

	/* This step accounts for any metric discontinuities at block boundaries
	   between grids at the same level */
	for(int level = fclaw_opt->maxlevel; level >= minlevel; level--)
	{
		/* Correct metric mismatches between patches of the same size */
		copy_at_blockbdry(glob,level,
						  read_parallel_patches,ghost_mode);

		time_sync_reset(glob, level, FCLAW2D_TIME_SYNC_RESET_SAMESIZE);

		/* Correct coarser grids with corrections from coarse/fine grids */
		if (level < fclaw_opt->maxlevel)
		{
			/* 'level' is the coarser level to be corrected */
			fine2coarse(glob,level,
			            read_parallel_patches,ghost_mode);

			/* Clear registers at coarse/fine interface */
			time_sync_reset(glob, level, FCLAW2D_TIME_SYNC_RESET_F2C);
		}
	}
	for(int level = fclaw_opt->maxlevel; level >= minlevel; level--)
	{
		/* Clear registers at physical level;  These are not used
     	   in sychronization, but will accumulate throughout simulation.  */
		time_sync_reset(glob, level, FCLAW2D_TIME_SYNC_RESET_PHYS);		
	}
}


void fclaw2d_time_sync(fclaw2d_global_t *glob, int minlevel, int maxlevel)
{


	fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_TIMESYNC]);

#if 0
	int time_interp = 0;

	/* --------------------------------------------------------------
		Send and receive patches.  This will carry out the exchange  
		for registers needed for conservation.
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



