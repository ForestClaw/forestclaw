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

#include <fclaw2d_global.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_partition.h>
#include <fclaw2d_exchange.h>
#include <fclaw2d_physical_bc.h>
#include <fclaw2d_regrid.h>
#include <fclaw2d_clawpatch.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* -----------------------------------------------------------------
   Initial grid
   ----------------------------------------------------------------- */
static
void cb_initialize (fclaw2d_domain_t *domain,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    void *user)
{
    fclaw2d_vtable_t vt;
    fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;

    fclaw2d_patch_data_new(domain,this_patch);
    fclaw2d_clawpatch_build(domain,this_patch,
                            this_block_idx,
                            this_patch_idx,
                            (void*) &build_mode);

    vt = fclaw2d_get_vtable(domain);
    vt.patch_initialize(domain,this_patch,this_block_idx,this_patch_idx);
}



/* -----------------------------------------------------------------
   Public interface
   ----------------------------------------------------------------- */

void fclaw2d_initialize (fclaw2d_domain_t **domain)
{
    int i;
    double t;

    char basename[BUFSIZ];
    const fclaw2d_vtable_t vt = fclaw2d_get_vtable(*domain);
    const amr_options_t *gparms = get_domain_parms(*domain);

    fclaw2d_domain_data_t* ddata = fclaw2d_domain_get_data(*domain);

    /* This mapping context is needed by fortran mapping functions */
    fclaw2d_map_context_t *cont = fclaw2d_domain_get_map_context(*domain);
    SET_CONTEXT(&cont);

    int minlevel = gparms->minlevel;
    int maxlevel = gparms->maxlevel;

    // This is where the timing starts.
    ddata->is_latest_domain = 1;
    for (i = 0; i < FCLAW2D_TIMER_COUNT; ++i) {
        fclaw2d_timer_init (&ddata->timers[i]);
    }

    /* set specific refinement strategy */
    fclaw2d_domain_set_refinement
      (*domain, gparms->smooth_refine, gparms->smooth_refine_level,
       gparms->coarsen_delay);

    /* start timing */
    fclaw2d_domain_barrier (*domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_INIT]);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_WALLTIME]);

    if (vt.problem_setup != NULL)
    {
        vt.problem_setup(*domain);
    }

    /* Get an initial domain */
    fclaw2d_domain_setup(NULL,*domain);

    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_BUILDREGRID]);
    fclaw2d_domain_iterate_level(*domain, minlevel, cb_initialize,
                                 (void *) NULL);
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_BUILDREGRID]);

    fclaw_bool time_interp = fclaw_false;
    t = fclaw2d_domain_get_time(*domain);
    fclaw2d_set_physical_bc(*domain,minlevel,t,time_interp);

    // VTK output during amrinit
    if (gparms->vtkout & 1) {
        // into timer
        fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_INIT]);
        fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);

        // output
        snprintf (basename, BUFSIZ, "%s_init_level_%02d",
                  gparms->prefix, minlevel);
        fclaw2d_output_write_vtk (*domain, basename);

        // out of timer
        fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);
        fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_INIT]);
    }

    // Refine as needed, one level at a time.
    if (minlevel < maxlevel)
    {
        int domain_init = 1;
        for (int level = minlevel; level < maxlevel; level++)
        {
            fclaw2d_domain_iterate_level(*domain, level,
                                         cb_fclaw2d_regrid_tag4refinement,
                                         (void *) &domain_init);

            // Construct new domain based on tagged patches.
            fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);
            int have_new_refinement = new_domain != NULL;

            // Domain data may go out of scope now.
            ddata = NULL;

            if (have_new_refinement)
            {
                fclaw2d_domain_setup(*domain,new_domain);
                ddata = fclaw2d_domain_get_data(new_domain);

                // Re-initialize new grids
                fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_BUILDREGRID]);
                fclaw2d_domain_iterate_adapted(*domain, new_domain,
                                               cb_fclaw2d_regrid_repopulate,
                                               (void *) &domain_init);
                fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_BUILDREGRID]);

                // Set boundary need ghost cell values so they are available
                // for using at tagging criteria, if necessary.
                int new_level = level+1;
                time_interp = fclaw_false;
                t = fclaw2d_domain_get_time(new_domain);
                fclaw2d_set_physical_bc(new_domain,new_level,t,time_interp);

                // free all memory associated with old domain
                fclaw2d_domain_reset(domain);
                *domain = new_domain;
                new_domain = NULL;

                // VTK output during amrinit
                if (gparms->vtkout & 1) {
                    // into timer
                    ddata = fclaw2d_domain_get_data (*domain);
                    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_INIT]);
                    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);

                    // output
                    snprintf (basename, BUFSIZ, "%s_init_level_%02d_adapt",
                              gparms->prefix, level);
                    fclaw2d_output_write_vtk (*domain, basename);

                    /* out of timer */
                    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);
                    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_INIT]);
                    ddata = NULL;
                }

                /* Repartition domain to new processors. Do I need this here?*/
                fclaw2d_partition_domain(domain, level);

                /* Set up ghost patches */
                fclaw2d_exchange_setup(*domain);
            }
            else
            {
                /* minlevel == maxlevel;  no refining necessary.  We have an initial
                   partition, so we don't need to partition a new domaoin. */

                fclaw2d_exchange_setup(*domain);
                break;
            }
        }
    }
    /* Print global minimum and maximum levels */
    fclaw_global_infof("Global minlevel %d maxlevel %d\n",
                (*domain)->global_minlevel, (*domain)->global_maxlevel);

    /* Stop timer */
    ddata = fclaw2d_domain_get_data(*domain);
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_INIT]);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
