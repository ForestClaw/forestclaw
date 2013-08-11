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



#include "amr_includes.H"

#include "amr_forestclaw.H"
#include "amr_mol.H"
#include "fclaw2d_solvers.H"


/* ----------------------------------------------------------
   Manage subcyling process
   ---------------------------------------------------------- */

static
double update_level_solution(fclaw2d_domain_t *domain,
                             int a_level,
                             fclaw2d_level_time_data *time_data)
{
    /* ToDo : Do we really need to pass in the entire time_data structure?
       Maybe this can be simplified considerably.
    */
    double t = time_data->t_level;
    double dt = time_data->dt;
    double cfl;

    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);

    /* Idea here is that the user may want to apply a single step routine,
       an MOL routine, or possibly both. */

    if (sf->use_single_step_update)
    {
        cfl = (sf->f_level_single_step)(domain,a_level,t,dt);
    }

    /* We may actually do both.  Just need to be sure that coarser level has taken a
       time step though */
    if (sf->use_mol_update)
    {
        cfl = fclaw2d_level_mol_step(domain,a_level,time_data,
                                     sf->f_level_ode_solver);
    }
    /* This is awkward. I should just return the maxcfl. */
    time_data->maxcfl = max(time_data->maxcfl,cfl);

    /* Okay, I am returning maxcfl.  But still need to clean up the time_data
       stuff.
    */
    return time_data->maxcfl;
}

static
double advance_level(fclaw2d_domain_t *domain,
                     const int a_level,
                     const int a_curr_fine_step,
                     double maxcfl,
                     subcycle_manager* a_time_stepper)
{
    // const amr_options_t *gparms = get_domain_parms(domain);
    fclaw_bool verbose = (fclaw_bool) a_time_stepper->verbosity();
    double t_level = a_time_stepper->level_time(a_level);

    int this_level = a_level;
    int coarser_level = a_level - 1;

    if (verbose)
    {
        cout << endl;
        cout << "Advancing level " << a_level << " from step " <<
            a_curr_fine_step << " at time " << t_level << endl;
    }

    /* -- Coming into this routine, all ghost cell information
          needed for an update of this level has been done.  So we can update
          immediately.

       -- All ghost cell exchanges needed to update coarser levels, whose
          last updated step is also 'a_curr_fine_step' have been done,
          and so coarser grids can all also be updated
    */

    fclaw2d_level_time_data_t time_data;

    time_data.t_level = t_level;
    time_data.t_initial = a_time_stepper->initial_time();
    time_data.dt = a_time_stepper->dt(a_level);
    time_data.fixed_dt = a_time_stepper->nosubcycle();

    if (verbose)
    {
        // cout << "Taking step on level " << a_level << " at time " << t_level << endl;
    }

    time_data.maxcfl = maxcfl;

#if 0
    /* Set some extra things needed by a multi-stage or implicit scheme. */
    time_data.is_coarsest = a_time_stepper->is_coarsest(a_level);
    if (!a_time_stepper->is_coarsest(a_level))
    {
        time_data.dt_coarse = a_time_stepper->dt(a_level-1);
    }
#endif

    /* ------------------------------------------------------------
       Advance this level from
       'a_curr_fine_step' to 'a_curr_fine_step + dt_level'
       ------------------------------------------------------------ */
    double cfl_step = update_level_solution(domain,this_level,&time_data);
    maxcfl = max(maxcfl,cfl_step);

    a_time_stepper->increment_step_counter(this_level);
    a_time_stepper->increment_time(this_level);


    /* Now advance coarser levels recursively.   If we are in the no-subcycle
       case, we will a take a time step of dt_fine (i.e. dt_level, where
       'level' is our current fine grid level).  If we are subcycling,
       then each time step will be the step size appropriate for that level.
       In this case, we are anticipating needing to do time interpolation to
       get ghost cells, and so as soon as we are finished with a coarse step,
       we will construct the time interpolated data.
    */

    if (!a_time_stepper->is_coarsest(this_level))
    {
        /* Advance coarser level, but only if coarse level and this level are time
           synchronized.  */
        int last_coarse_step = a_time_stepper->last_step(coarser_level);
        if (last_coarse_step == a_curr_fine_step)
        {
            double cfl_step = advance_level(domain,coarser_level,last_coarse_step,
                                            maxcfl,a_time_stepper);
            maxcfl = max(maxcfl,cfl_step);

            if (!a_time_stepper->nosubcycle())
            {
                /* Time interpolate this data for a future exchange with finer grid */
                int coarse_inc = a_time_stepper->step_inc(coarser_level);
                int new_curr_step = a_time_stepper->last_step(this_level);
                double alpha = double(new_curr_step % coarse_inc)/coarse_inc;
                if (verbose)
                {
                    cout << "Time interpolating level " << coarser_level <<
                        " using alpha = " << alpha << endl;
                }

                timeinterp(domain,coarser_level,alpha);
            }
        }
    }


    if (verbose)
    {
        cout << "Advance on level " << a_level << " done at time " <<
            a_time_stepper->level_time(a_level) << endl << endl;
    }

    return time_data.maxcfl;  // Maximum from level iteration
}


/* ------------------------------------------------------------
   This does a complete exchange over all levels.

   -- Parallel ghost patches are exchange at all levels
   -- Every level exchanges ghost cells with other patches
      at that level
   -- Every finer level exchanges with a coarser level
   -- No time interpolation is assumed, as all levels are time
      synchronized at this point.
   -- This the only routine that is called for the non-subcycled
      case.
   -- All levels will be updated in next update step, regardless of
      whether we are in the subcycled or non-subcycled case.

  The reason for two separate ghost cell exchange routines is that
  the logic here is considerably simpler than for the partial
  update used in intermediate steps in the subcycled case.
  ------------------------------------------------------------- */
void update_ghost_all_levels(fclaw2d_domain_t* domain,
                             fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);

    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw_bool verbose = gparms->verbosity;

    if (verbose)
    {
        cout << "Exchanging ghost patches across all levels " << endl;
    }

    exchange_ghost_patch_data_all(domain);

    /* Each level should exchange with other patches at that level */
    int minlevel = domain->global_minlevel;
    int maxlevel = domain->global_maxlevel;

    fclaw_bool time_interp = fclaw_false;
    double t = get_domain_time(domain);

    for(int level = maxlevel; level >= minlevel; level--)
    {
        level_exchange(domain,level,FCLAW2D_TIMER_EXCHANGE);
    }

    /* Each level should exchange with coarser patches */
    for(int level = maxlevel; level > minlevel; level--)
    {
        exchange_with_coarse(domain,level,t,time_interp,
                             FCLAW2D_TIMER_EXCHANGE);
    }

    for(int level = maxlevel; level >= minlevel; level--)
    {
        set_phys_bc(domain,level,t,time_interp);
    }

    // Stop timing
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_start (&ddata->timers[running]);
    }
}



/* ----------------------------------------------------------
   This does the intermediate ghost patch and ghost cell
   exchange needed by either the subcycled case or the non
   subcycled case.

   Assumptions going into this routine :
   -- coarse_level > global_minlevel
   -- Data at coarse_level-1 is time interpolated data
   -- This routine is never called in the non-subcycled case
   -- This should be more 'lightweight' than doing a complete
      exchange (across all levels) at each fine grid time
      step.  But it should be equivalent.

   ---------------------------------------------------------- */
static
void update_ghost_partial(fclaw2d_domain_t* domain, int coarse_level,
                          int fine_level, subcycle_manager *a_time_stepper,
                          fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);

    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw_bool verbose = gparms->verbosity;

    int time_interp_level = coarse_level - 1;

    if (verbose)
    {
        cout << "Exchanging ghost patches from levels " <<
            coarse_level << " to " << fine_level << endl;
        if (!a_time_stepper->nosubcycle())
        {
            cout << "Time interpolated level is " <<
                time_interp_level << endl;
        }
    }

    /* Make available patches from levels coarse to fine */
    set_boundary_patch_ptrs(domain,time_interp_level, fine_level);

    /* Do parallel ghost exchange */
    exchange_ghost_patch_data_levels(domain,time_interp_level,fine_level);

    /* -------------------------------------------------------
       Do ghost cell exchange.
       ------------------------------------------------------- */
    double t = get_domain_time(domain); /* needed for phys. bc */
    // fclaw_bool time_interp = fclaw_false;
    /* Do all level exchanges first */
    for(int level = fine_level; level >= coarse_level; level--)
    {
        level_exchange(domain,level,FCLAW2D_TIMER_EXCHANGE);
    }


    /* Then do exchanges with coarser level */
    fclaw_bool time_interp_vec[fine_level-coarse_level+1];
    for(int level = fine_level; level > coarse_level; level--)
    {
        time_interp_vec[level] = fclaw_false;
    }

    /* coarse level exchanges with a time interpolated level
       time_interp_level = coarse_level-1
    */
    time_interp_vec[coarse_level] = fclaw_true;
    for(int level = fine_level; level >= coarse_level; level--)
    {
        exchange_with_coarse(domain,level,t,time_interp_vec[level],
                             FCLAW2D_TIMER_EXCHANGE);
    }


    fclaw_bool time_interp = fclaw_false;
    for (int level = fine_level; level >= coarse_level; level--)
    {
        set_phys_bc(domain,level,t,time_interp);
    }

    // Stop timing
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_start (&ddata->timers[running]);
    }
}



/* -------------------------------------------------------------
   Main routine : Called from amrrun.cpp
   ------------------------------------------------------------- */

double advance_all_levels(fclaw2d_domain_t *domain,
                          subcycle_manager *a_time_stepper)
{
    // Start timer for advancing the levels
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
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
        double cfl_step = advance_level(domain,maxlevel,nf,maxcfl,a_time_stepper);
        maxcfl = max(cfl_step,maxcfl);
        int last_step = a_time_stepper->last_step(maxlevel);
        if (!a_time_stepper->nosubcycle() && last_step < n_fine_steps)
        {
            /* Find time interpolated level and do ghost patch exchange
             and ghost cell exchange for next update.
            */
            int time_interp_level = maxlevel-1;
            while (last_step % a_time_stepper->step_inc(time_interp_level) == 0)
            {
                time_interp_level--;
            }
            /* This exchange includes an exchange between level maxlevel-n+1
               and the time interpolated patch at level maxlevel-n.
            */
            update_ghost_partial(domain,time_interp_level+1,maxlevel,
                                 a_time_stepper,FCLAW2D_TIMER_ADVANCE);
        }
    }
    /* Do a complete update.  This is needed even if we are regridding
       before the  next step.  In this case, the ghost cell values can
       be used in tagging for refining.  If we regrid, this will have to
       be called a second time.
       Idea : If we know we are regridding before the next step, we could
       skip this update, and make it clear to the user that ghost cell
       values are not available for determining refinement critera.
    */
    update_ghost_all_levels(domain,FCLAW2D_TIMER_ADVANCE);

    // Stop the timer
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_ADVANCE]);
    ++ddata->count_amr_advance;

    return maxcfl;
}
