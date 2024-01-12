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

#include <fclaw2d_regrid.h>

#include <fclaw2d_options.h>

#include <fclaw2d_convenience.h>   /* p4est domain, patch handling routines */
#include <forestclaw2d.h>          /* Needed for patch_relation_t data */

#include <fclaw_gauges.h>

#include <fclaw2d_global.h>
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_partition.h>
#include <fclaw2d_exchange.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>


/* This is also called from fclaw2d_initialize, so is not made static */
void cb_fclaw2d_regrid_tag4refinement(fclaw2d_domain_t *domain,
									  fclaw2d_patch_t *this_patch,
									  int this_block_idx,
									  int this_patch_idx,
									  void *user)
{
    int refine_patch, maxlevel, level;
    const fclaw_options_t* fclaw_opt;

    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;
    int domain_init = *((int*) g->user);

    fclaw_opt = fclaw2d_get_options(g->glob);

    maxlevel = fclaw_opt->maxlevel;
    level = this_patch->level;

    if (level < maxlevel)
    {
        refine_patch  =
            fclaw2d_patch_tag4refinement(g->glob,this_patch,this_block_idx,
                                         this_patch_idx, domain_init);
        if (refine_patch == 1)
        {
            fclaw2d_patch_mark_refine(domain, this_block_idx, this_patch_idx);
        }
    }
}

/* Tag family for coarsening */

void cb_regrid_tag4coarsening(fclaw2d_domain_t *domain,
							  fclaw2d_patch_t *fine_patches,
							  int blockno, int fine0_patchno,
							  void *user)
{
    const fclaw_options_t* fclaw_opt;
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;

    fclaw_opt = fclaw2d_get_options(g->glob);
    int domain_init = *((int*) g->user);

    int minlevel = fclaw_opt->minlevel;

    int level = fine_patches[0].level;

    if (level > minlevel)
    {
        int family_coarsened = 1;
        family_coarsened = fclaw2d_patch_tag4coarsening(g->glob,&fine_patches[0],
                                                        blockno, fine0_patchno,
                                                        domain_init);
        if (family_coarsened == 1)
        {
            int igrid;
            for (igrid = 0; igrid < 4; igrid++)
            {
                int fine_patchno = fine0_patchno + igrid;
                fclaw2d_patch_mark_coarsen(domain,blockno, fine_patchno);
            }
        }
    }
}


/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */

void cb_fclaw2d_regrid_repopulate(fclaw2d_domain_t * old_domain,
								  fclaw2d_patch_t * old_patch,
								  fclaw2d_domain_t * new_domain,
								  fclaw2d_patch_t * new_patch,
								  fclaw2d_patch_relation_t newsize,
								  int blockno,
								  int old_patchno,
								  int new_patchno,
								  void *user)
{
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;
    int domain_init = *((int*) g->user);

    fclaw2d_domain_data_t *ddata_old = fclaw2d_domain_get_data (old_domain);
    fclaw2d_domain_data_t *ddata_new = fclaw2d_domain_get_data (new_domain);
    fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;

    if (newsize == FCLAW2D_PATCH_SAMESIZE)
    {
        FCLAW_ASSERT(0 <= blockno && blockno < new_domain->num_blocks);
        FCLAW_ASSERT(0 <= new_patchno && new_patchno < new_domain->local_num_patches);
        new_patch->user = old_patch->user;
        old_patch->user = NULL;
        ++ddata_old->count_delete_patch;
        ++ddata_new->count_set_patch;
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        fclaw2d_patch_t *fine_siblings = new_patch;
        fclaw2d_patch_t *coarse_patch = old_patch;

        int i;
        for (i = 0; i < 4; i++)
        {
            fclaw2d_patch_t *fine_patch = &fine_siblings[i];
            int fine_patchno = new_patchno + i;
            /* Reason for the following two lines: the glob contains the old domain which is incremented in ddata_old 
               but we really want to increment the new domain. This will be fixed! */
            --ddata_old->count_set_patch;
            ++ddata_new->count_set_patch;

            fclaw2d_patch_build(g->glob,fine_patch,blockno,
                                fine_patchno,(void*) &build_mode);
            if (domain_init)
            {
                fclaw2d_patch_initialize(g->glob,fine_patch,blockno,fine_patchno);//new_domain
            }
        }

        if (!domain_init)
        {
            int coarse_patchno = old_patchno;
            int fine_patchno = new_patchno;

            fclaw2d_patch_interpolate2fine(g->glob,coarse_patch,fine_siblings,
                                           blockno,coarse_patchno,fine_patchno);//new_domain
        }
        /* used to pass in old_domain */
        fclaw2d_patch_data_delete(g->glob,coarse_patch);
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {

#if 0      
        if (domain_init)
        {
            /* We now do coarsening at the initial refinement */
            fclaw_debugf("fclaw2d_regrid.cpp (repopulate): We shouldn't end up here\n");
            exit(0);
        }
#endif        

        /* Old grids are the finer grids;  new grid is the coarsened grid */
        fclaw2d_patch_t *fine_siblings = old_patch;
        int fine_patchno = old_patchno;

        fclaw2d_patch_t *coarse_patch = new_patch;
        int coarse_patchno = new_patchno;
        
        /* Reason for the following two lines: the glob contains the old domain which is incremented in ddata_old 
           but we really want to increment the new domain. This will be fixed! */
        --ddata_old->count_set_patch;
        ++ddata_new->count_set_patch;
        
        if (domain_init)
        {
            fclaw2d_patch_build(g->glob,coarse_patch,blockno,
                                coarse_patchno,(void*) &build_mode);
            fclaw2d_patch_initialize(g->glob,coarse_patch,blockno,coarse_patchno);
        }
        else
        {
            /* Area (and possibly other things) should be averaged to coarse grid. */
            fclaw2d_patch_build_from_fine(g->glob,fine_siblings,coarse_patch,
                                          blockno,coarse_patchno,fine_patchno,
                                          build_mode);
            /* Average the solution. Does this need to be customizable? */
            fclaw2d_patch_average2coarse(g->glob,fine_siblings,coarse_patch,
                                        blockno,fine_patchno,coarse_patchno);

        }
        int i;
        for(i = 0; i < 4; i++)
        {
            fclaw2d_patch_t* fine_patch = &fine_siblings[i];
            /* used to pass in old_domain */
            fclaw2d_patch_data_delete(g->glob,fine_patch);
        }
    }
    else
    {
        fclaw_global_essentialf("cb_adapt_domain : newsize not recognized\n");
        exit(1);
    }
    fclaw2d_patch_neighbors_reset(new_patch);
}

/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */
void fclaw2d_regrid(fclaw2d_global_t *glob)
{
    fclaw2d_domain_t** domain = &glob->domain;
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_REGRID]);

    fclaw_global_infof("Regridding domain\n");

    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_REGRID_TAGGING]);
    /* First determine which families should be coarsened. */
    int domain_init = 0;
    fclaw2d_global_iterate_families(glob, cb_regrid_tag4coarsening,
                                    (void *) &domain_init);

    fclaw2d_global_iterate_patches(glob, cb_fclaw2d_regrid_tag4refinement,
                                   (void *) &domain_init);

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_REGRID_TAGGING]);

    /* Rebuild domain if necessary */
    /* Will return be NULL if no refining was done */

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_REGRID]);
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_ADAPT_COMM]);
    fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);

    int have_new_refinement = new_domain != NULL;

    if (have_new_refinement)
    {
        /* allocate memory for user patch data and user domain data in the new
           domain;  copy data from the old to new the domain. */
        fclaw2d_domain_setup(glob, new_domain);
    }

    /* Stop the new timer (copied from old timer) */
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_ADAPT_COMM]);
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_REGRID]);

    if (have_new_refinement)
    {
        fclaw_global_infof(" -- Have new refinement\n");

        /* Average to new coarse grids and interpolate to new fine grids */
        fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_REGRID_BUILD]);
        fclaw2d_global_iterate_adapted(glob, new_domain,
                                       cb_fclaw2d_regrid_repopulate,
                                       (void *) &domain_init);
        fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_REGRID_BUILD]);

        /* free memory associated with old domain */
        fclaw2d_domain_reset(glob);
        *domain = new_domain;
        new_domain = NULL;

        /* Repartition for load balancing.  Second arg (mode) for vtk output */
        fclaw2d_partition_domain(glob,FCLAW2D_TIMER_REGRID);

        /* Set up ghost patches. Communication happens for indirect ghost exchanges. */


        /* This includes timers for building patches and (exclusive) communication */
        fclaw2d_exchange_setup(glob,FCLAW2D_TIMER_REGRID);

        /* Get new neighbor information.  This is used to short circuit
           ghost filling procedures in some cases */
        fclaw2d_regrid_set_neighbor_types(glob);

        /* Update ghost cells.  This is needed because we have new coarse or fine
           patches without valid ghost cells.   Time_interp = 0, since we only
           only regrid when all levels are time synchronized. */
        int minlevel = (*domain)->global_minlevel;
        int maxlevel = (*domain)->global_maxlevel;
        int time_interp = 0;
        double sync_time = glob->curr_time;
        fclaw2d_ghost_update(glob,
                             minlevel,
                             maxlevel,
                             sync_time,
                             time_interp,
                             FCLAW2D_TIMER_REGRID);

        ++glob->count_amr_new_domain;
    }
    else
    {
        /* We updated all the ghost cells when leaving advance, so don't need to do
           it here */
    }

    /* User defined */
    fclaw2d_after_regrid(glob, have_new_refinement);

    /* Only if gauges count > 0 */
    fclaw_locate_gauges(glob);

    /* Stop timer.  Be sure to use timers from new grid, if one was
       created */
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_REGRID]);

    /* Count calls to this function */
    ++glob->count_amr_regrid;
}

static
void cb_set_neighbor_types(fclaw2d_domain_t *domain,
						   fclaw2d_patch_t *this_patch,
						   int blockno,
						   int patchno,
						   void *user)
{
	int iface, icorner;

	for (iface = 0; iface < 4; iface++)
	{
		int rproc[2];
		int rblockno;
		int rpatchno[2];
		int rfaceno;

		fclaw2d_patch_relation_t neighbor_type =
		fclaw2d_patch_face_neighbors(domain,
									 blockno,
									 patchno,
									 iface,
									 rproc,
									 &rblockno,
									 rpatchno,
									 &rfaceno);

		fclaw2d_patch_set_face_type(this_patch,iface,neighbor_type);
	}

	for (icorner = 0; icorner < 4; icorner++)
	{
		int rproc_corner;
		int cornerpatchno;
		int cornerblockno;
		int rcornerno;
		fclaw2d_patch_relation_t neighbor_type;

		int has_corner_neighbor =
		fclaw2d_patch_corner_neighbors(domain,
									   blockno,
									   patchno,
									   icorner,
									   &rproc_corner,
									   &cornerblockno,
									   &cornerpatchno,
									   &rcornerno,
									   &neighbor_type);

		fclaw2d_patch_set_corner_type(this_patch,icorner,neighbor_type);
		if (!has_corner_neighbor)
		{
			fclaw2d_patch_set_missing_corner(this_patch,icorner);
		}
	}
	fclaw2d_patch_neighbors_set(this_patch);
}

/* Set neighbor type : samesize, halfsize, or doublesize */
void fclaw2d_regrid_set_neighbor_types(fclaw2d_global_t *glob)
{
	fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_NEIGHBOR_SEARCH]);
	fclaw2d_global_iterate_patches(glob,cb_set_neighbor_types,NULL);
	fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_NEIGHBOR_SEARCH]);
}
