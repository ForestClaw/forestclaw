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

static void update_level_solution(fclaw2d_domain_t *domain,
                                  int a_level,
                                  fclaw2d_level_time_data *time_data)
{
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
    time_data->maxcfl = max(time_data->maxcfl,cfl);
}

static
double advance_level(fclaw2d_domain_t *domain,
                   const int& a_level,
                   const int& a_curr_fine_step,
                   subcycle_manager* a_time_stepper)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw_bool verbose = (fclaw_bool) a_time_stepper->verbosity();
    double t_level = a_time_stepper->level_time(a_level);

    double maxcfl = 0;
    if (verbose)
    {
        cout << endl;
        cout << "Advancing level " << a_level << " from step " <<
            a_curr_fine_step << " at time " << t_level << endl;
    }

    if (!a_time_stepper->can_advance(a_level,a_curr_fine_step))
    {
        if (!a_time_stepper->level_exchange_done(a_level))
        {
            /* Level exchange should have been done right after solution update. */
            printf("Error (advance_level) : Level exchange at level %d not done at time step %d\n",
                   a_level,a_curr_fine_step);
            exit(1);
        }
        if (!a_time_stepper->exchanged_with_coarser(a_level))
        {
            int last_coarse_step = a_time_stepper->last_step(a_level-1);
            if (verbose)
            {
                cout << " --> Exchange between fine level " << a_level <<
                    " and coarse level " << a_level-1 << " not done at step "
                     << a_curr_fine_step << endl;
            }
            if (a_curr_fine_step == last_coarse_step)
            {
                /* Levels are time synchronized and we can do a simple coarse/fine
                   exchange without time interpolation or first advancing the
                   coarser level */
                if (verbose)
                {
                    cout << " ----> Coarse and fine level are time synchronized; doing exchange"
                        " without time interpolation" << endl;
                }
                double alpha = 0;  // No time interpolation

                /* Boundary conditions on finer grid should be set so that they can
                   be used in setting corners on coarser grids, if needed. */
                set_phys_bc(domain,a_level,t_level);

                exchange_with_coarse(domain,a_level,t_level,alpha);

                // set_phys_bc(domain,a_level,t_level);
                a_time_stepper->increment_coarse_exchange_counter(a_level);
                a_time_stepper->increment_fine_exchange_counter(a_level-1);

                if (a_time_stepper->nosubcycle() && !a_time_stepper->is_coarsest(a_level))
                {
                    /* non-subcycled case : this advance is only needed because we
                       want to advance the coarser grids, not because we need time
                       interpolated boundary conditions. */
                    if (verbose)
                    {
                        cout << " ----> Non-subcycled case " << endl;
                        cout << " ----> Making recursive call to advance_level for level "
                             << a_level-1 << endl;
                    }
                    maxcfl = advance_level(domain,a_level-1,last_coarse_step,a_time_stepper);
                    if (verbose)
                    {
                        cout << " ----> Returning from recursive call at level "
                             << a_level << endl;
                    }
                }
            }
            else
            {
                /* We may need time interpolated boundary conditions. */
                if (verbose)
                {
                    cout << " --> Coarse and fine level are not time synchronized; doing exchange "
                        "with time interpolation" << endl;
                }
                if ((a_curr_fine_step > last_coarse_step))
                {
                    /* subcycled case : a_curr_fine_step will only be greater than
                       last_coarse_step if we haven't yet advanced the coarse grid to a time level
                       beyond the current fine grid level. */
                    if (verbose)
                    {
                        cout << " ----> Subcycled case " << endl;
                        cout << " ----> Making recursive call to advance_level for level "
                             << a_level-1 << endl;
                    }
                    maxcfl = advance_level(domain,a_level-1,last_coarse_step,a_time_stepper);
                    if (verbose)
                    {
                        cout << " ----> Returning from recursive call at level "
                             << a_level << endl;
                    }
                }
                else
                {
                    /* This will only happen if refratio > 2.  In that case, we have more than
                       one fine grid step which can use the same two coarse grid time levels
                       to get time interpolated boundary conditions. */
                }
                if (!a_time_stepper->nosubcycle())
                {
                    if (verbose)
                    {
                        cout << " --> Doing time interpolatation from coarse grid at level "
                             << a_level-1 << endl;
                    }
                    int refratio = gparms->refratio;

                    /* (1) a_curr_fine_step > last_coarse_step : we just advanced the coarse grid
                       and can now apply time interpolated boundary conditions, or

                       (2) a_curr_fine_step < last_coarse_step : we advanced the coarse
                       grid in a previous step but we still have to do time interpolation (this
                       can only happen if refratio > 2) */

                    double  alpha = double(a_curr_fine_step)/refratio;
                    if (verbose)
                    {
                        cout << " --> Doing time interpolatation from coarse grid at level "
                             << a_level-1 << " using alpha = " << alpha << endl;
                    }
                    exchange_with_coarse(domain,a_level,t_level,alpha);
                    set_phys_bc(domain,a_level,t_level);
                    a_time_stepper->increment_coarse_exchange_counter(a_level);
                }
            }
        }
    }
    fclaw2d_level_time_data_t time_data;

    time_data.t_level = t_level;
    time_data.t_initial = a_time_stepper->initial_time();
    time_data.dt = a_time_stepper->dt(a_level);
    time_data.fixed_dt = a_time_stepper->nosubcycle();

    if (verbose)
    {
        cout << "Taking step on level " << a_level << " at time " << t_level << endl;
    }

    /* Need to figure out what to do if our time step doesn't depend on
       a cfl condition */
    time_data.maxcfl = maxcfl;

    /* Set some extra things needed by a multi-stage or implicit scheme. */
    time_data.is_coarsest = a_time_stepper->is_coarsest(a_level);
    if (!a_time_stepper->is_coarsest(a_level))
    {
        time_data.dt_coarse = a_time_stepper->dt(a_level-1);
    }

    /* ------------------------------------------------------------
       Advance this level from
       'a_curr_fine_step' to 'a_curr_fine_step + dt_level'
       ------------------------------------------------------------ */
    update_level_solution(domain,a_level,&time_data);

    a_time_stepper->increment_step_counter(a_level);
    a_time_stepper->increment_time(a_level);

    level_exchange(domain,a_level);
    set_phys_bc(domain,a_level,t_level);

    set_phys_bc(domain,a_level,t_level);
    a_time_stepper->increment_level_exchange_counter(a_level);

    if (verbose)
    {
        cout << "Advance on level " << a_level << " done at time " <<
            a_time_stepper->level_time(a_level) << endl << endl;
    }

    return time_data.maxcfl;  // Maximum from level iteration
}

/* -------------------------------------------------------------
   Main routine : Called from amrrun.cpp
   ------------------------------------------------------------- */

double advance_all_levels(fclaw2d_domain_t *domain,
                        subcycle_manager *a_time_stepper)
{
    // 'n_fine_steps' is the number of steps we must take on the finest level to equal one
    // step on the coarsest non-empty level, i.e. minlevel.
    int minlevel = a_time_stepper->minlevel();
    int n_fine_steps = a_time_stepper->step_inc(minlevel); // equal 1 in the non-subcycled case.
    int maxlevel = a_time_stepper->maxlevel();
    double maxcfl = 0;
    for(int nf = 0; nf < n_fine_steps; nf++)
    {
        double cfl_step = advance_level(domain,maxlevel,nf,a_time_stepper);
        maxcfl = max(cfl_step,maxcfl);
    }
    return maxcfl;
}
