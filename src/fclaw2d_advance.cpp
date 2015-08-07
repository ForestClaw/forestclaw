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


#include <fclaw2d_forestclaw.h>
#include <fclaw_base.h>
#include <fclaw2d_advance.hpp>
#include <fclaw2d_timeinterp.h>
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_update_single_step.h>

#include <math.h>



/* ----------------------------------------------------------
   Manage subcyling process
   ---------------------------------------------------------- */

static
double update_level_solution(fclaw2d_domain_t *domain,
                             int a_level,
                             fclaw2d_level_time_data *time_data)
{
    /* ToDo : Do we really need to pass in the entire time_data
       structure? Maybe this can be simplified considerably.
    */
    double t = time_data->t_level;
    double dt = time_data->dt;
    double cfl;

    /* Note : there might not be any grids at this level */
    cfl = fclaw2d_update_single_step(domain,a_level,t,dt);

    /* This needs to be cleaned up a bit */
    time_data->maxcfl = fmax(time_data->maxcfl,cfl);

    return cfl;
}

static
double advance_level(fclaw2d_domain_t *domain,
                     const int a_level,
                     const int a_curr_fine_step,
                     double maxcfl,
                     subcycle_manager* a_time_stepper)
{
    double t_level = a_time_stepper->level_time(a_level);

    int this_level = a_level;
    int coarser_level = a_level - 1;

    fclaw_global_infof("Advancing level %d from step %d at time %12.6e\n",
                       this_level,a_curr_fine_step,t_level);


    /* -- Coming into this routine, all ghost cell information
          needed for an update of this level has been done.  So we can
          update immediately.

       -- All ghost cell exchanges needed to update coarser levels,
          whose last updated step is also 'a_curr_fine_step' have been
          done, and so coarser grids can all also be updated
    */

    /* The use of time_data here could be cleaned up a bit, but I am
       leaving everything for now, in case I later decide to do more
       with the MOL approach.
    */
    fclaw2d_level_time_data_t time_data;

    time_data.t_level = t_level;
    time_data.t_initial = a_time_stepper->initial_time();
    time_data.dt = a_time_stepper->dt(this_level);
    time_data.fixed_dt = a_time_stepper->nosubcycle();
    time_data.maxcfl = maxcfl;

    /* ------------------------------------------------------------
       Advance this level from 'a_curr_fine_step' to
       'a_curr_fine_step + dt_level'  This returns 0 if there are
       no grids at this level.
       ------------------------------------------------------------ */
    double cfl_step = update_level_solution(domain,this_level,
                                            &time_data);
    maxcfl = fmax(maxcfl,cfl_step);

    fclaw_global_infof("------ Max CFL on level %d is %12.4e " \
                       " (using dt = %12.4e)\n",this_level,cfl_step,time_data.dt);

    a_time_stepper->increment_step_counter(this_level);
    a_time_stepper->increment_time(this_level);


    /* Advance coarser levels recursively.   If we are in the no-subcycle
       case, we will a take a time step of dt_fine (i.e. dt_level, where
       'level' is our current fine grid level).  If we are subcycling,
       then each time step will be the step size appropriate for that
       level. In this case, we are anticipating needing to do time
       interpolation to get ghost cells, and so as soon as we are
       finished with a coarse step, we will construct the time
       interpolated data.
    */

    if (!a_time_stepper->is_coarsest(this_level))
    {
        /* Advance coarser level, but only if coarse level and this
           level are time synchronized.  */
        int last_coarse_step = a_time_stepper->last_step(coarser_level);
        if (last_coarse_step == a_curr_fine_step)
        {
            double cfl_step = advance_level(domain,coarser_level,
                                            last_coarse_step,
                                            maxcfl,a_time_stepper);
            maxcfl = fmax(maxcfl,cfl_step);

            if (!a_time_stepper->nosubcycle())
            {
                /* Time interpolate this data for a future exchange with finer grid */
                int coarse_inc =
                    a_time_stepper->step_inc(coarser_level);
                int new_curr_step =
                    a_time_stepper->last_step(this_level);
                double alpha =
                    double(new_curr_step % coarse_inc)/coarse_inc;

                fclaw_global_infof("Time interpolating level %d using alpha = %5.2f\n",
                                   coarser_level,alpha);

                /* This happens even if there are no grids at the
                   time interp level */
                fclaw2d_timeinterp(domain,coarser_level,alpha);
            }
        }
    }

    fclaw_global_infof("Advance on level %d done at time %12.6e\n\n",
                       a_level,a_time_stepper->level_time(a_level));

    return maxcfl;  // Maximum from level iteration
}

/* -------------------------------------------------------------
   Main routine : Called from fclaw2d_run.
   ------------------------------------------------------------- */

double advance_all_levels(fclaw2d_domain_t *domain,
                          subcycle_manager *a_time_stepper)
{
    // Start timer for advancing the levels
    fclaw2d_domain_data_t* ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_ADVANCE]);

    /* These are global minimum and maximum values */
    int minlevel = a_time_stepper->minlevel();
    int maxlevel = a_time_stepper->maxlevel();

    /* Number of steps we need to take on the finest level to equal one
       step on the coarsest level.  In the non-subcycled case, this
       'n_fine_steps' is equal to 1.  In either case, we have
       n_fine_steps = 2^(maxlevel-minlevel)
    */
    int n_fine_steps = a_time_stepper->step_inc(minlevel);

    /* Keep track of largest cfl over all grid updates */
    double maxcfl = 0;
    for(int nf = 0; nf < n_fine_steps; nf++)
    {
        double cfl_step = advance_level(domain,maxlevel,nf,maxcfl,
                                        a_time_stepper);
        maxcfl = fmax(cfl_step,maxcfl);
        int last_step = a_time_stepper->last_step(maxlevel);
        if (!a_time_stepper->nosubcycle() && last_step < n_fine_steps)
        {
            /* Find time interpolated level and do ghost patch exchange
               and ghost cell exchange for next update. */
            int time_interp_level = maxlevel-1;
            while (last_step %
                   a_time_stepper->step_inc(time_interp_level) == 0)
            {
                time_interp_level--;
            }

            /* This exchange includes an exchange between level
               maxlevel-n+1 and the time interpolated patch at level
               maxlevel-n. */

            /* coarsest level is a time interpolated level */
            int time_interp = 1;
#if 0
            fclaw2d_ghost_update(domain,minlevel,maxlevel,
                                 time_interp,FCLAW2D_TIMER_ADVANCE);
#endif
            double sync_time = a_time_stepper->level_time(last_step);

            fclaw2d_ghost_update(domain,
                                 time_interp_level+1,
                                 maxlevel,
                                 sync_time,
                                 time_interp,
                                 FCLAW2D_TIMER_ADVANCE);
        }
    }
    /* Do a complete update.  This is needed even if we are regridding
       before the  next step.  In this case, the ghost cell values can
       be used in tagging for refining.  If we regrid, this will have to
       be called a second time.
       Idea : If we know we are regridding before the next step, we could
       skip this update, and make it clear to the user that ghost cell
       values are not available for determining refinement criteria.
    */
    fclaw_global_infof("Advance is done with coarse grid step at " \
                      " time %12.6e\n",a_time_stepper->initial_time());

    /* We have to do this here, since we may need ghost cells for interpolating
       coarse grids to newly created fine grids.  If we have a new refinement,
       we follow the regridding step with a second update_ghost */

    int time_interp = 0;

    /* Physical time needed for physical boundary conditions */
    double sync_time =  a_time_stepper->level_time(n_fine_steps);
    fclaw2d_ghost_update(domain,minlevel,maxlevel,sync_time,
                         time_interp,FCLAW2D_TIMER_ADVANCE);

    // Stop the timer
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_ADVANCE]);
    ++ddata->count_amr_advance;

    return maxcfl;
}
