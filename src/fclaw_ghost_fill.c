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

/*
 * IDEA 1: add boolean about on_parallel_boundary vs. interior.
 *         1. Do the parallel boundary work here.
 *            To the least possible amount of work before sending.
 *             * copy_samelevel
 *             * average_ghost_f2c
 *             * fill_physical_ghost
 *               (which subset of cells shall do it?)
 *               (maybe we can skip doing it before sending altogether?)
 *            More specifically, only copy/average across faces
 *            between two parallel boundary patches, not corners.
 *            Argue in paper why no corners needed for 5-point stencils.
 *         2. Start sending: ghost_exchange_begin.
 *         3. Do all the work on interior patches and whatever has not
 *            been computed between parallel boundary patches.
 *            Do the most possible amount of work while communicating.
 *             * copy_samelevel
 *             * average_ghost_f2c
 *             * fill_physical_ghost
 *               (how much of this can we possibly to here: maximize this.)
 *             * interpolate_ghost_c2f for interior patches and
 *               faces between parallel boundary and interior patches.
 *         4. Recieve messages: ghost_exchange_end.
 *         5. Work on receive buffers / parallel patches with remote data.
 *             * copy_samelevel
 *             * average_ghost_f2c
 *             * fill_physical_ghost
 *             * interpolate_ghost_c2f
 *         6. Repeat 5. for whatever parts of patches that are not done yet.
 *            Integrate this with 5. so we don't loop over patches twice.
 *
 * All of this goes into the paper as detailed algorithm with explanations.
 */


/** \file
 * Fill ghost cells.
 *
 *
 **/


#include <fclaw_ghost_fill.h>
#include <fclaw_corner_neighbors.h>
#include <fclaw_edge_neighbors.h>
#include <fclaw_face_neighbors.h>

#include <fclaw_global.h>
#include <fclaw_timer.h>
#include <fclaw_patch.h>
#include <fclaw_partition.h>
#include <fclaw_exchange.h>
#include <fclaw_domain.h>
#include <fclaw_physical_bc.h>


/* -------------------------------------------------
   Basic routines - operate on a single level
   ----------------------------------------------- */

//static fclaw2d_patch_iterator_t patch_iterator;

typedef struct fclaw2d_ghost_fill_wrap_info
{
	fclaw_ghost_fill_parallel_mode_t ghost_mode;
	fclaw_patch_callback_t cb_fill;
	void *user;
} fclaw2d_ghost_fill_wrap_info_t;

static void cb_parallel_wrap(fclaw_domain_t* domain,
							 fclaw_patch_t* this_patch,
							 int this_block_idx,
							 int this_patch_idx,
							 void *user);


static void cb_interface_wrap(fclaw_domain_t* domain,
							  fclaw_patch_t* this_patch,
							  int this_block_idx,
							  int this_patch_idx,
							  void *user)
{
	if (fclaw_patch_on_coarsefine_interface(this_patch))
	{
		cb_parallel_wrap(domain,
						 this_patch,
						 this_block_idx,
						 this_patch_idx,
						 user);
	}
}

/* these routines filter patches according to whether we are working on
   patches on the parallel boundary only, or on interior patches.
*/
/* Passing in a function pointer is somehow slower ... */
static void cb_parallel_wrap(fclaw_domain_t* domain,
							 fclaw_patch_t* this_patch,
							 int this_block_idx,
							 int this_patch_idx,
							 void *user)
{
	fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;  
	
	fclaw2d_ghost_fill_wrap_info_t w = *((fclaw2d_ghost_fill_wrap_info_t*) s->user);
	int on_boundary = fclaw_patch_on_parallel_boundary(this_patch);
	if (((w.ghost_mode == FCLAW_BOUNDARY_GHOST_ONLY) && on_boundary) ||
		((w.ghost_mode == FCLAW_BOUNDARY_INTERIOR_ONLY) && !on_boundary) ||
		w.ghost_mode == FCLAW_BOUNDARY_ALL)
	{
		fclaw_global_iterate_t g;
		g.glob = s->glob;
		g.user = w.user;    
		
		w.cb_fill(domain,this_patch,this_block_idx,this_patch_idx,&g);
	}
}


static
void copy2ghost(fclaw_global_t *glob,
				int level,
				int time_interp,
				int read_parallel_patches,
				fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
    struct fclaw2d_ghost_fill_wrap_info parallel_mode;
    fclaw_exchange_info_t e_info;
    e_info.exchange_type = FCLAW_COPY;
    e_info.grid_type = FCLAW_IS_COARSE;
    e_info.time_interp = time_interp;
    e_info.read_parallel_patches = read_parallel_patches;

    parallel_mode.ghost_mode = ghost_mode;
    parallel_mode.user = (void*) &e_info;

    parallel_mode.cb_fill = fclaw_face_fill_cb;
    fclaw_global_iterate_level(glob, level, cb_parallel_wrap,
                         (void *) &parallel_mode);

    if(glob->domain->refine_dim == 3)
    {
        /* edge exchanges */
        parallel_mode.cb_fill = fclaw_edge_fill_cb;
        fclaw_global_iterate_level(glob, level, cb_parallel_wrap,
                              (void *) &parallel_mode);
    }

    /* corner exchanges */
    parallel_mode.cb_fill = fclaw_corner_fill_cb;
    fclaw_global_iterate_level(glob, level, cb_parallel_wrap,
                          (void *) &parallel_mode);
}


static
void average2ghost(fclaw_global_t *glob,
				   int coarse_level,
				   int time_interp,
				   int read_parallel_patches,
				   fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw2d_ghost_fill_wrap_info_t parallel_mode;

	fclaw_exchange_info_t e_info;
	e_info.time_interp = time_interp; /* Does this matter here? */
	e_info.read_parallel_patches = read_parallel_patches;
	e_info.exchange_type = FCLAW_AVERAGE;
	e_info.has_fine_grid_neighbor = 0;

	/* Only update ghost cells at local boundaries */
	e_info.grid_type = FCLAW_IS_COARSE;

	parallel_mode.user = (void*) &e_info;
	parallel_mode.ghost_mode = ghost_mode;

    parallel_mode.cb_fill = fclaw_face_fill_cb;
    fclaw_global_iterate_level(glob, coarse_level,
                   cb_interface_wrap, (void *) &parallel_mode);

	if(glob->domain->refine_dim ==3)
    {
        /* Edge average */
        parallel_mode.cb_fill = fclaw_edge_fill_cb;
        fclaw_global_iterate_level(glob, coarse_level, cb_interface_wrap,
                       (void *) &parallel_mode);
    }

    /* Corner average */
    parallel_mode.cb_fill = fclaw_corner_fill_cb;
    fclaw_global_iterate_level(glob, coarse_level, cb_interface_wrap,
                   (void *) &parallel_mode);

	if (read_parallel_patches)
	{
		/* Second pass : average from local fine grids to remote coarse grids. These
		   coarse grids might be needed for interpolation later. */
		e_info.grid_type = FCLAW_IS_FINE;

		int fine_level = coarse_level + 1;

		/* Face average */
		parallel_mode.cb_fill = fclaw_face_fill_cb;
		fclaw_global_iterate_level(glob, fine_level,
									 cb_interface_wrap, (void *) &parallel_mode);

		if(glob->domain->refine_dim == 3)
	    {
	        /* Interpolate to edges at parallel boundaries from coarse grid
	           ghost patches */
	        parallel_mode.cb_fill = fclaw_edge_fill_cb;
	        fclaw_global_iterate_level(glob, fine_level, cb_interface_wrap,
	                                     (void *) &parallel_mode);
	    }

		parallel_mode.cb_fill = fclaw_corner_fill_cb;
		fclaw_global_iterate_level(glob, fine_level, cb_interface_wrap,
									 (void *) &parallel_mode);
	}
}

static
void interpolate2ghost(fclaw_global_t *glob,
					   int coarse_level,
					   int time_interp,
					   int read_parallel_patches,
					   fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	struct fclaw2d_ghost_fill_wrap_info parallel_mode;
	fclaw_exchange_info_t e_info;
	e_info.time_interp = time_interp;
#if 0
	e_info.level = coarse_level;
#endif
    e_info.level = -1;
    e_info.exchange_type = FCLAW_INTERPOLATE;

    /* ----------------------------------------------------------
       First pass - look for fine grids to interpolate to. This
       should include include the time interpolated level.
       ---------------------------------------------------------- */

    e_info.grid_type = FCLAW_IS_COARSE;
    e_info.read_parallel_patches = read_parallel_patches;

    parallel_mode.ghost_mode = ghost_mode;
    parallel_mode.user = (void*) &e_info;

    /* Face interpolate */
    parallel_mode.cb_fill = fclaw_face_fill_cb;
    fclaw_global_iterate_level(glob,coarse_level, cb_interface_wrap,
                                         (void *) &parallel_mode);

    if(glob->domain->refine_dim ==3)
    {
        /* Edge interpolate */
        parallel_mode.cb_fill = fclaw_edge_fill_cb;
        fclaw_global_iterate_level(glob, coarse_level, cb_interface_wrap,
                                  (void *) &parallel_mode);
    }

    /* Corner interpolate */
    parallel_mode.cb_fill = fclaw_corner_fill_cb;
    fclaw_global_iterate_level(glob,coarse_level, cb_interface_wrap,
                  (void *) &parallel_mode);

    if (read_parallel_patches)
    {
        /* -----------------------------------------------------
           Second pass - Iterate over local fine grids, looking
           for remote coarse grids we can use to fill in BCs at
           fine grid ghost cells along the parallel boundary
           ----------------------------------------------------- */

        e_info.grid_type = FCLAW_IS_FINE;

        /* Interpolate to faces at parallel boundaries from coarse grid ghost
           patches */
        int fine_level = coarse_level + 1;

        parallel_mode.cb_fill = fclaw_face_fill_cb;
        fclaw_global_iterate_level(glob, fine_level, cb_interface_wrap,
                                   (void *) &parallel_mode);

		if(glob->domain->refine_dim == 3)
	    {
	        parallel_mode.cb_fill = fclaw_edge_fill_cb;
	        fclaw_global_iterate_level(glob, fine_level, cb_interface_wrap,
	                                     (void *) &parallel_mode);
	    }

        /* Interpolate to corners at parallel boundaries from coarse grid
           ghost patches */
        parallel_mode.cb_fill = fclaw_corner_fill_cb;
        fclaw_global_iterate_level(glob, fine_level, cb_interface_wrap,
                                   (void *) &parallel_mode);
    }
}


static
void setphysical(fclaw_global_t *glob,
				 int level,
				 double sync_time,
				 int time_interp,
				 fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	struct fclaw2d_ghost_fill_wrap_info parallel_mode;

	fclaw_physical_time_info_t t_info;
	t_info.level_time = sync_time;
	t_info.time_interp = time_interp;

	parallel_mode.ghost_mode = ghost_mode;
	parallel_mode.user = (void*) &t_info;

	parallel_mode.cb_fill = cb_fclaw_physical_set_bc;
	fclaw_global_iterate_level(glob, level, cb_parallel_wrap,
								 (void *) &parallel_mode);
}


/* -------------------------------------------------
   Loop over all levels
   ----------------------------------------------- */

#if 0
static
void time_sync(fclaw2d_global_t* glob,
			   int minlevel,
			   int maxlevel,
			   int time_interp,
			   int read_parallel_patches,
			   fclaw2d_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_GHOSTFILL_COPY]);
	int level;

	/* Copy between grids that are at the same level. */
	for(level = maxlevel; level >= minlevel; level--)
	{
		int time_interp = 0;
		time_sync_copy(glob,level,time_interp,
				   read_parallel_patches,ghost_mode);
	}
	for(level = maxlevel; level >= minlevel; level--)
	{
		int time_interp = 0;
		time_sync_fine_to_coarse(glob,level,time_interp,
								read_parallel_patches,ghost_mode);
	}
	fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_GHOSTFILL_COPY]);
}
#endif


static
void copy_samelevel(fclaw_global_t* glob,
					int minlevel,
					int maxlevel,
					int time_interp,
					int read_parallel_patches,
					fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_COPY]);
	int level;

	/* Copy between grids that are at the same level. */
	for(level = maxlevel; level >= minlevel; level--)
	{
		int time_interp = 0;
		copy2ghost(glob,level,time_interp,
				   read_parallel_patches,ghost_mode);
	}
	if (time_interp)
	{
		int time_interp_level = minlevel-1;
		copy2ghost(glob,time_interp_level,time_interp,
				   read_parallel_patches,ghost_mode);
	}
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_COPY]);
}


/*
 * Fill in coarse grid ghost cells by averaging or copying  from neighboring fine grids.
 */
static
void average_fine2coarse_ghost(fclaw_global_t *glob,
							   int mincoarse,
							   int maxcoarse,
							   int time_interp,
							   int read_parallel_patches,
							   fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_AVERAGE]);

	int level;

	/* Average fine grids to coarse grid ghost cells */
	for(level = maxcoarse; level >= mincoarse; level--)
	{
		/* Don't do any time interpolation yet */
		int time_interp = 0;
		average2ghost(glob,level,time_interp,
					  read_parallel_patches,ghost_mode);
	}
	if (time_interp)
	{
		/* Average fine grid to coarse time interpolated level.  Time interpolated
		 faces need correct values.  */
		int time_interp_level = mincoarse-1;
		average2ghost(glob,time_interp_level,time_interp,
					  read_parallel_patches,ghost_mode);
	}
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_AVERAGE]);
}

static
void interpolate_coarse2fine_ghost(fclaw_global_t* glob,
								   int mincoarse,
								   int maxcoarse,
								   int time_interp,
								   int read_parallel_patches,
								   fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_INTERP]);

	int level;

	/* Interpolate from coarse grid to fine grid ghost */
	for(level = maxcoarse; level >= mincoarse; level--)
	{
		/* No need to interpolate to coarse time-interpolated grids */
		int time_interp = 0;
		interpolate2ghost(glob,level,time_interp,read_parallel_patches,ghost_mode);
	}
	if (time_interp)
	{
		/* interpolate from time interpolated level */
		int time_interp_level = mincoarse-1;
		interpolate2ghost(glob,time_interp_level,time_interp,
						  read_parallel_patches,ghost_mode);
	}
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_INTERP]);
}


static
void fill_physical_ghost(fclaw_global_t* glob,
						 int minlevel,
						 int maxlevel,
						 double sync_time,
						 int time_interp,
						 fclaw_ghost_fill_parallel_mode_t ghost_mode)
{
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_PHYSBC]);

	int level;

	for(level = maxlevel; level >= minlevel; level--)
	{
		int time_interp = 0;
		setphysical(glob,level,sync_time,time_interp,ghost_mode);
	}
	if (time_interp)
	{
		int time_interp_level = minlevel-1;
		setphysical(glob,
					time_interp_level,
					sync_time,
					time_interp,
					ghost_mode);
	}
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_PHYSBC]);
}

/* -----------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------------*/

void fclaw_ghost_update_nonasync(fclaw_global_t* glob,
								   int minlevel,
								   int maxlevel,
								   double sync_time,
								   int time_interp,
								   fclaw_timer_names_t running)
{
	if (running != FCLAW_TIMER_NONE) {
		fclaw_timer_stop (&glob->timers[running]);
	}
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL]);

	fclaw_global_infof("Exchanging ghost patches across all levels\n");

#if 0
	/* uncomment this if debugging ghost cell interpolation */
	fclaw_global_essentialf("WARNING : compute_slopes set to average\n");
#endif

	/* ---------------------------------------------------------
	   Get coarse grid ghost cells ready to use for interpolation.
	   Coarse ghost cells on ghost patches are not updated in
	   this step.  Ghost patches do not yet have valid data, and
	   so boundary patches will have to be updated after the exchange.
	   ---------------------------------------------------------- */

	if (time_interp)
	{
		int time_interp_level = minlevel - 1;
		fclaw_global_infof("Time interpolated level is %d\n",   \
						   time_interp_level);
	}

	/* If minlevel == maxlevel, then maxcoarse < mincoase. In this
	   case, loops involving averaging and interpolation will be
	   skipped and we only copy between patches on
	   same levels. */

	int mincoarse = minlevel;
	int maxcoarse = maxlevel-1;   /* maxlevel >= minlevel */

    /* --------------------------------------------------------------
    Do work we have do before sending
    ------------------------------------------------------------*/
    /* Copy/average ghost cells in local patches at the parallel boundary.
        This is needed so that when these boundary patches get sent to other
        processors as ghost patches, they have valid ghost cells if needed
        for interpolation.*/
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP1]);

    fclaw_ghost_fill_parallel_mode_t parallel_mode 
                        = FCLAW_BOUNDARY_GHOST_ONLY;
    int read_parallel_patches = 0;


    /* Average */
    average_fine2coarse_ghost(glob,mincoarse,maxcoarse,
                              time_interp,
                              read_parallel_patches,
                              parallel_mode);

    /* Copy */
    copy_samelevel(glob,minlevel,maxlevel,time_interp,
                   read_parallel_patches,parallel_mode);


    /* This is needed when the parallel boundary intersects the physical
       boundary.  In this case, a coarse grid ghost patch might
       need to have physical boundary conditions in order to interpolate
       to a fine grid local patch. Note that we don't need to fill in
       ghost cells on the finest level grids, since these will never be
        used for interpolation. */
    fill_physical_ghost(glob,
                        mincoarse,
                        maxcoarse,
                        sync_time,
                        time_interp,
                        parallel_mode);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP1]);

	/* --------------------------------------------------------------
		Send and receive ghost patches
	------------------------------------------------------------*/
	fclaw_exchange_ghost_patches_begin(glob,minlevel,maxlevel,time_interp,
										 FCLAW_TIMER_GHOSTFILL);

	fclaw_exchange_ghost_patches_end(glob,minlevel,maxlevel,time_interp,
									   FCLAW_TIMER_GHOSTFILL);

	/* -------------------------------------------------------------
		Indirect neighbor exchange
	------------------------------------------------------------- */
	fclaw_face_neighbor_ghost(glob,minlevel,maxlevel,time_interp);
	fclaw_corner_neighbor_indirect(glob,minlevel,maxlevel, time_interp);


	/* -------------------------------------------------------------
		Do everything in one set of sequential steps
	------------------------------------------------------------- */

	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP3]);

    read_parallel_patches = 1;

    parallel_mode =
           FCLAW_BOUNDARY_ALL;    

	/* Average */
	average_fine2coarse_ghost(glob,mincoarse,maxcoarse, time_interp,
							  read_parallel_patches,parallel_mode);

	/* Copy */
	copy_samelevel(glob,minlevel,maxlevel,time_interp,
				   read_parallel_patches,parallel_mode);

	/* Physical */
	fill_physical_ghost(glob,
						minlevel,
						maxlevel,
						sync_time,
						time_interp,
						parallel_mode);

	/* Interpolate */
	interpolate_coarse2fine_ghost(glob,mincoarse, maxcoarse,
								  time_interp,read_parallel_patches,
								  parallel_mode);
	/* Physical */
	fill_physical_ghost(glob,
						minlevel,
						maxlevel,
						sync_time,
						time_interp,
						FCLAW_BOUNDARY_ALL);
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP3]);

	// Stop timing
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL]);
	if (running != FCLAW_TIMER_NONE)
	{
		fclaw_timer_start (&glob->timers[running]);
	}
}


void fclaw_ghost_update_async(fclaw_global_t* glob,
								int minlevel,
								int maxlevel,
								double sync_time,
								int time_interp,
								fclaw_timer_names_t running)
{
	if (running != FCLAW_TIMER_NONE) {
		fclaw_timer_stop (&glob->timers[running]);
	}
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL]);

	fclaw_global_infof("Exchanging ghost patches across all levels\n");

#if 0
	/* uncomment this if debugging ghost cell interpolation */
	fclaw_global_essentialf("WARNING : compute_slopes set to average\n");
#endif

	/* ---------------------------------------------------------
	   Get coarse grid ghost cells ready to use for interpolation.
	   Coarse ghost cells on ghost patches are not updated in
	   this step.  Ghost patches do not yet have valid data, and
	   so boundary patches will have to be updated after the exchange.
	   ---------------------------------------------------------- */

	if (time_interp)
	{
		int time_interp_level = minlevel - 1;
		fclaw_global_infof("Time interpolated level is %d\n",   \
						   time_interp_level);
	}

	/* If minlevel == maxlevel, then maxcoarse < mincoase. In this
	   case, loops involving averaging and interpolation will be
	   skipped and we only copy between patches on
	   same levels. */

	int mincoarse = minlevel;
	int maxcoarse = maxlevel-1;   /* maxlevel >= minlevel */

	fclaw_ghost_fill_parallel_mode_t parallel_mode =
		   FCLAW_BOUNDARY_ALL;
	int read_parallel_patches;

	/* --------------------------------------------------------------
	Do work we have do before sending
	------------------------------------------------------------*/
	/* Copy/average ghost cells in local patches at the parallel boundary.
		This is needed so that when these boundary patches get sent to other
		processors as ghost patches, they have valid ghost cells if needed
		for interpolation.*/
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP1]);

	parallel_mode = FCLAW_BOUNDARY_GHOST_ONLY;
	read_parallel_patches = 0;


	/* Average */
	average_fine2coarse_ghost(glob,mincoarse,maxcoarse,
							  time_interp,
							  read_parallel_patches,
							  parallel_mode);

	/* Copy */
	copy_samelevel(glob,minlevel,maxlevel,time_interp,
				   read_parallel_patches,parallel_mode);


	/* This is needed when the parallel boundary intersects the physical
	   boundary.  In this case, a coarse grid ghost patch might
	   need to have physical boundary conditions in order to interpolate
	   to a fine grid local patch. Note that we don't need to fill in
	   ghost cells on the finest level grids, since these will never be
		used for interpolation. */
	fill_physical_ghost(glob,
						mincoarse,
						maxcoarse,
						sync_time,
						time_interp,
						parallel_mode);

	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP1]);

	/* --------------------------------------------------------------
		Start send ...
	------------------------------------------------------------*/
	fclaw_exchange_ghost_patches_begin(glob,minlevel,maxlevel,time_interp,
										 FCLAW_TIMER_GHOSTFILL);

	/* --------------------------------------------------------------
		Finish exchanges in the interior of the grid.
	------------------------------------------------------------*/

	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP2]);
	parallel_mode = FCLAW_BOUNDARY_INTERIOR_ONLY;

	/* Average */
	average_fine2coarse_ghost(glob,mincoarse,maxcoarse,
							  time_interp,
							  read_parallel_patches,
							  parallel_mode);

	/* Copy */
	copy_samelevel(glob,minlevel,maxlevel,time_interp,
				   read_parallel_patches,parallel_mode);

	/* Physical ghost */
	fill_physical_ghost(glob,
						mincoarse,
						maxcoarse,
						sync_time,
						time_interp,
						parallel_mode);

	/* Fine grids that are adjacent to boundary patches don't get
		ghost regions filled in that overlap boundary patch */

	/* Interpolate */
	interpolate_coarse2fine_ghost(glob,mincoarse, maxcoarse,
								  time_interp,
								  read_parallel_patches,
								  parallel_mode);

	/* minfine to maxfine?  */
	int minfine = mincoarse+1;
	int maxfine = maxlevel;

	/* Physical ghost */
	fill_physical_ghost(glob,
						minfine,
						maxfine,
						sync_time,
						time_interp,
						parallel_mode);

	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP2]);

	/* -------------------------------------------------------------
		Receive ghost patches ...
	------------------------------------------------------------- */

	fclaw_exchange_ghost_patches_end(glob,minlevel,maxlevel,time_interp,
									   FCLAW_TIMER_GHOSTFILL);

	/* -------------------------------------------------------------
		Loop over ghost patches to find indirect neighbors and do
		any necessary face exchanges.

		Note : There is no special timer for this call, but checks
		show that ghostfill-(step1+step2+step3+comm) << 1
	------------------------------------------------------------- */
	fclaw_face_neighbor_ghost(glob,minlevel,maxlevel,time_interp);
	fclaw_corner_neighbor_indirect(glob,minlevel,maxlevel, time_interp);

	/* -------------------------------------------------------------
		Repeat above, but now with parallel ghost cells.
	------------------------------------------------------------- */

	fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP3]);
	parallel_mode = FCLAW_BOUNDARY_GHOST_ONLY;
	read_parallel_patches = 1;

	/* Average */
	average_fine2coarse_ghost(glob,mincoarse,maxcoarse, time_interp,
							  read_parallel_patches,parallel_mode);

	/* Copy */
	copy_samelevel(glob,minlevel,maxlevel,time_interp,
				   read_parallel_patches,parallel_mode);

	/* Physical */
	fill_physical_ghost(glob,
						minlevel,
						maxlevel,
						sync_time,
						time_interp,
						parallel_mode);

	/* Interpolate */
	interpolate_coarse2fine_ghost(glob,mincoarse, maxcoarse,
								  time_interp,read_parallel_patches,
								  parallel_mode);

	/*  This needs to be over all patches.  The situation that can arise is
		that a coarse grid boundary patch adjacent to a fine grid non-boundary
		patch. Both have a physical boundary.  In this case, the following steps
		happens at the coarse/fine interface :

		(1) Ghost only (before exchange) : Fine grid averages to coarse grid;
		Physical BC applied to coarse grid.

		(2) Interior only : Fine grid is not picked up in an interpolation step,
			since in interpolation, we only sweep over coarse grids, looking for fine
			grid neighbors.  Applying physical BC to this fine grid leaves invalid values
			in fine grid ghost cells at boundary.

		(3) Ghost only (after exchange) : Coarse grid boundary patch interpolates
		   to adjacent fine grid.  But now, we need to apply physical BCs to fine grid.
		   So we sweep over all grids to apply phys. BCs.  This is overkill, since
		   only those fine grids with neighboring coarse grid patches on the parallel boundary
		   are affected, but it is hard to see how to avoid this without some tedious
		   checking.
	*/

	/* Physical */
	fill_physical_ghost(glob,
						minlevel,
						maxlevel,
						sync_time,
						time_interp,
						FCLAW_BOUNDARY_ALL);
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL_STEP3]);

	// Stop timing
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTFILL]);
	if (running != FCLAW_TIMER_NONE)
	{
		fclaw_timer_start (&glob->timers[running]);
	}
}



/* -----------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------------*/

void fclaw_ghost_update(fclaw_global_t* glob,
						  int minlevel,
						  int maxlevel,
						  double sync_time,
						  int time_interp,
						  fclaw_timer_names_t running)
{
    
	int async = 1;
	if (async == 0)
	{
		fclaw_ghost_update_nonasync(glob,minlevel,maxlevel,sync_time,time_interp,running);
	}
	else
	{
		/* About 10% faster? */
		fclaw_ghost_update_async(glob,minlevel,maxlevel,sync_time,time_interp,running);
	}
}


