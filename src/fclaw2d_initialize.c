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

#include <fclaw2d_forestclaw.h>

#include <fclaw_global.h>
#include <fclaw_options.h>

#include <fclaw2d_convenience.h>

#include <fclaw_gauges.h>

#include <fclaw2d_partition.h>
#include <fclaw2d_exchange.h>
#include <fclaw2d_physical_bc.h>
#include <fclaw2d_regrid.h>
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_map.h>
#include <fclaw_patch.h>
#include <fclaw_domain.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

/* -----------------------------------------------------------------
   Initial grid
   ----------------------------------------------------------------- */
static
void cb_initialize (fclaw_domain_t *domain,
					fclaw_patch_t *this_patch,
					int this_block_idx,
					int this_patch_idx,
					void *user)
{
	fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;

	fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;

	fclaw_patch_build(g->glob,this_patch,
						this_block_idx,
						this_patch_idx,
						&build_mode);
	fclaw_patch_initialize(g->glob,this_patch,this_block_idx,this_patch_idx);
}



/* -----------------------------------------------------------------
   Public interface
   ----------------------------------------------------------------- */
void fclaw2d_initialize(fclaw_global_t *glob)
{
	fclaw_domain_t** domain = &glob->domain;

    int time_interp = 0;
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

	/* This mapping context is needed by fortran mapping functions */
	fclaw2d_map_context_t *cont = fclaw_global_get_map_2d(glob);
	FCLAW_MAP_SET_CONTEXT(&cont);

	int maxthreads = 0;

#if defined(_OPENMP)
	maxthreads = omp_get_max_threads();
#endif
	
    fclaw_global_essentialf("Max threads set to %d\n",maxthreads);

    int minlevel = fclaw_opt->minlevel;
    int maxlevel = fclaw_opt->maxlevel;

    /* Initialize all timers */
    int i;
    for (i = 0; i < FCLAW_TIMER_COUNT; ++i) {
        fclaw_timer_init (&glob->timers[i]);
    }

    /* start timing */
    fclaw2d_domain_barrier (*domain);
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_WALLTIME]);
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_INIT]);

    /* User defined problem setup */
    fclaw2d_problem_setup(glob);

    /* set specific refinement strategy */
    fclaw2d_domain_set_refinement
        (*domain, fclaw_opt->smooth_refine, fclaw_opt->smooth_level,
         fclaw_opt->coarsen_delay);


    /* ------------------------------------------------
       Set up initial domain.

       This needs to be set as if it were going to be used
       for updating.
       ------------------------------------------------ */

    /* Get an initial domain */
    fclaw_domain_setup(glob,*domain);

    /* Initialize patches on uniformly refined level minlevel */
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_REGRID_BUILD]);
    fclaw_global_iterate_level(glob, minlevel, cb_initialize,
                                 NULL);
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_REGRID_BUILD]);

    /* Set up ghost patches */
    fclaw2d_exchange_setup(glob,FCLAW_TIMER_INIT);

    /* This is normally called from regrid */
    fclaw2d_regrid_set_neighbor_types(glob);

    /* We need a user option here to set ghost values after initialization */
    if (fclaw_opt->init_ghostcell)
    {
        fclaw2d_ghost_update(glob,(*domain)->global_minlevel,
                             (*domain)->global_maxlevel,0.0,
                             time_interp,FCLAW_TIMER_INIT);
    }

#if 0
    fclaw2d_physical_set_bc(glob,(*domain)->global_minlevel,
                            0.0,fclaw_opt->dt_initial,time_interp);
#endif                            

    /* ------------------------------------------------
       Done with initial setup.
       ------------------------------------------------ */


    /* ------------------------------------------------
       Build up an initial refinement.
       ------------------------------------------------ */
    if (minlevel < maxlevel)
    {
        int domain_init = 1;
        int level;
        for (level = minlevel; level <= maxlevel; level++)
        {
            fclaw_timer_start (&glob->timers[FCLAW_TIMER_REGRID_TAGGING]);

            fclaw_global_iterate_families(glob, cb_regrid_tag4coarsening,
                                            (void*) &domain_init);

            fclaw_global_iterate_patches(glob,cb_fclaw2d_regrid_tag4refinement,
                                         &domain_init);
            
            fclaw_timer_stop (&glob->timers[FCLAW_TIMER_REGRID_TAGGING]);

            // Construct new domain based on tagged patches.
            fclaw_timer_stop (&glob->timers[FCLAW_TIMER_INIT]);
            fclaw_timer_start (&glob->timers[FCLAW_TIMER_ADAPT_COMM]);
            fclaw_domain_t *new_domain = fclaw2d_domain_adapt(*domain);

            int have_new_refinement = new_domain != NULL;

            if (have_new_refinement)
            {
                /* Have to get a new ddata */
                fclaw_domain_setup(glob,new_domain);
            }

            fclaw_timer_stop (&glob->timers[FCLAW_TIMER_ADAPT_COMM]);
            fclaw_timer_start (&glob->timers[FCLAW_TIMER_INIT]);

            if (have_new_refinement)
            {
                fclaw_global_infof(" -- Have new initial refinement\n");

                /* Re-initialize new grids.   Ghost cell values needed for
                   interpolation have already been set by initialization */
                fclaw_timer_start (&glob->timers[FCLAW_TIMER_REGRID_BUILD]);
                fclaw_global_iterate_adapted(glob, new_domain,
                                               cb_fclaw2d_regrid_repopulate,
                                               (void *) &domain_init);

                fclaw_timer_stop (&glob->timers[FCLAW_TIMER_REGRID_BUILD]);

                /* free all memory associated with old domain */
                fclaw_domain_reset(glob);
                *domain = new_domain;
                new_domain = NULL;

                /* Repartition domain to new processors.    */
                fclaw2d_partition_domain(glob,FCLAW_TIMER_INIT);

                /* Set up ghost patches.  This probably doesn't need to be done
                   each time we add a new level. */
                fclaw2d_exchange_setup(glob,FCLAW_TIMER_INIT);
                
                /* This is normally called from regrid, once the initial domain
                   has been set up */
                fclaw2d_regrid_set_neighbor_types(glob);
            }
            else
            {
                fclaw_global_infof(" -- No new initial refinement\n");
                /* We don't have a new refinement, and so can break out of level loop */
                break;
            }
        }  /* Level loop (minlevel --> maxlevel) */
    }

    if (fclaw_opt->init_ghostcell)
    {
        fclaw2d_ghost_update(glob,(*domain)->global_minlevel,
                             (*domain)->global_maxlevel,0.0,
                             time_interp,FCLAW_TIMER_INIT);
    }

    fclaw2d_diagnostics_initialize(glob);
    fclaw_locate_gauges(glob);

    fclaw_after_regrid(glob);

    /* Print global minimum and maximum levels */
    fclaw_global_infof("Global minlevel %d maxlevel %d\n",
                (*domain)->global_minlevel, (*domain)->global_maxlevel);

    /* Stop timer */
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_INIT]);
}
