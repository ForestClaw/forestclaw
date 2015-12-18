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

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_update_single_step.h>


#include <fclaw2d_advance.hpp>

#if 0
#include <fclaw2d_timestep_counters.hpp>
#endif


/* ----------------------------------------------------------
   Manage subcyling process
   ---------------------------------------------------------- */

static
double update_level_solution(fclaw2d_domain_t *domain,
                             int level,
                             double t, double dt)
{
    /* There might not be any grids at this level */
    double cfl = fclaw2d_update_single_step(domain,level,t,dt);

    return cfl;
}

static
double advance_level(fclaw2d_domain_t *domain,
                     const int level,
                     const int curr_fine_step,
                     double maxcfl,
                     fclaw2d_timestep_counters* ts_counter)
{
    const amr_options_t* gparms = get_domain_parms(domain);

    double t_level = ts_counter->level_time(level);
    double dt_level = ts_counter->dt(level);

    int this_level = level;
    int coarser_level = level - 1;

    fclaw_global_infof("Advancing level %d from step %d at time %12.6e\n",
                       this_level,curr_fine_step,t_level);

    double cfl_step = update_level_solution(domain,this_level,t_level,dt_level);

    maxcfl = fmax(maxcfl,cfl_step);

    fclaw_global_infof("------ Max CFL on level %d is %12.4e " \
                       " (using dt = %12.4e)\n",this_level,cfl_step,dt_level);

    ts_counter->increment_step_counter(this_level);
    ts_counter->increment_time(this_level);

    if (this_level > domain->local_minlevel)
    {
        int last_coarse_step = ts_counter->last_step(coarser_level);
        if (last_coarse_step == curr_fine_step)
        {
            double cfl_step = advance_level(domain,coarser_level,
                                            last_coarse_step,
                                            maxcfl,ts_counter);
            maxcfl = fmax(maxcfl,cfl_step);

            if (gparms->subcycle)
            {
                double alpha = ts_counter->compute_alpha(this_level);

                fclaw_global_infof("Time interpolating level %d using alpha = %5.2f\n",
                                   coarser_level,alpha);

                fclaw2d_timeinterp(domain,coarser_level,alpha);
            }
        }
    }

    fclaw_global_infof("Advance on level %d done at time %12.6e\n\n",
                       level,ts_counter->level_time(level));

    return maxcfl;  /* Maximum from level iteration */
}

/* -------------------------------------------------------------
   Main routine : Called from fclaw2d_run.
   ------------------------------------------------------------- */

double fclaw2d_advance_all_levels(fclaw2d_domain_t *domain,
                                  double t_curr, double dt)
{
    fclaw2d_domain_data_t* ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_ADVANCE]);

    const amr_options_t* gparms = get_domain_parms(domain);

    /* Take steps on the finest level grids present anywhere */
    fclaw2d_timestep_counters ts_counter;
    ts_counter.define(domain,gparms,t_curr,dt);

    /* Advance all grids that are present somewhere (on any proc) */
    int minlevel = domain->global_minlevel;
    int maxlevel = domain->global_maxlevel;

    /* Keep track of largest cfl over all grid updates */
    double maxcfl = 0;

    int n_fine_steps = ts_counter.steps(maxlevel);
    for(int nf = 0; nf < n_fine_steps; nf++)
    {
        /* Coarser levels get updated recursively */
        double cfl_step = advance_level(domain,maxlevel,nf,maxcfl,&ts_counter);

        maxcfl = fmax(cfl_step,maxcfl);
        int last_step = nf+1;
        if (last_step < n_fine_steps)
        {
            /* Do intermediate ghost cell exchanges */
            double sync_time = ts_counter.level_time(maxlevel);
            if (gparms->subcycle)
            {
                /* Find time interpolated level and do ghost patch exchange
                   and ghost cell exchange for next update. */
                int time_interp_level = ts_counter.timeinterp_level(maxlevel);
                int time_interp = 1;
                fclaw2d_ghost_update(domain,
                                     time_interp_level+1,
                                     maxlevel,
                                     sync_time,
                                     time_interp,
                                     FCLAW2D_TIMER_ADVANCE);
            }
            else
            {
                int time_interp = 0;
                fclaw2d_ghost_update(domain,
                                     minlevel,
                                     maxlevel,
                                     sync_time,
                                     time_interp,
                                     FCLAW2D_TIMER_ADVANCE);
            }
        }
    }
    fclaw_global_infof("Advance is done with coarse grid step at " \
                      " time %12.6e\n",t_curr);

    double sync_time =  ts_counter.level_time(maxlevel);
    FCLAW_ASSERT(fabs(sync_time - (t_curr+dt)) < 1e-13);
    int time_interp = 0;
    fclaw2d_ghost_update(domain,minlevel,maxlevel,sync_time,
                         time_interp,FCLAW2D_TIMER_ADVANCE);

    // Stop the timer
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_ADVANCE]);
    ++ddata->count_amr_advance;

    return maxcfl;
}
