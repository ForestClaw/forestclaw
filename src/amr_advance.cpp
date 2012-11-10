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

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "fclaw2d_convenience.h"
#include "fclaw_defs.H"

class ClawPatch;

static
void cb_advance_patch(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    const amr_options_t* gparms = ddata->amropts;

    ClawPatch *cp = get_clawpatch(this_patch);
    fclaw2d_level_time_data_t *time_data = (fclaw2d_level_time_data_t *) user;

    double dt = time_data->dt;
    double t = time_data->t;

    int level = this_patch->level;

    double maxcfl_grid = cp->step_noqad(t,dt,level,*gparms);
    time_data->maxcfl = max(maxcfl_grid,time_data->maxcfl);
}

static
double advance_level(fclaw2d_domain_t *domain,
                   const int& a_level,
                   const int& a_curr_fine_step,
                   subcycle_manager* a_time_stepper)
{
    bool verbose = (bool) a_time_stepper->m_verbosity;
    double t_level = a_time_stepper->current_time(a_level);

    double maxcfl = 0;
    if (verbose)
    {
        cout << endl;
        cout << "Advancing level " << a_level << " from step " <<
            a_curr_fine_step << endl;
    }

    if (!a_time_stepper->can_advance(a_level,a_curr_fine_step))
    {
        if (!a_time_stepper->level_exchange_done(a_level))
        {
            // Level exchange should have been done right after solution update.
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
                // Levels are time synchronized and we can do a simple coarse/fine
                // exchange without time interpolation or advancing the coarser level
                if (verbose)
                {
                    cout << " ----> Coarse and fine level are time synchronized; doing exchange"
                        " without time interpolation" << endl;
                }
                exchange_with_coarse(domain,a_level);

                set_phys_bc(domain,a_level,t_level);
                a_time_stepper->increment_coarse_exchange_counter(a_level);
                a_time_stepper->increment_fine_exchange_counter(a_level-1);

                if (a_time_stepper->nosubcycle() && !a_time_stepper->is_coarsest(a_level))
                {
                    // non-subcycled case : this advance is a bit gratuitous, because we
                    // don't need it to advance the fine grid;  rather, we put this advance
                    // here as a way to get advances on the coarser levels.
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
                if (verbose)
                {
                    cout << " --> Coarse and fine level are not time synchronized; doing exchange "
                        "with time interpolation" << endl;
                }
                if ((a_curr_fine_step > last_coarse_step))
                {
                    // subcycled case : a_curr_fine_step will only be greater than
                    // last_coarse_step if we haven't yet advanced the coarse grid to a time level
                    // beyond the current fine grid level.
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
                if (!a_time_stepper->nosubcycle())
                {
                    if (verbose)
                    {
                        cout << " --> Doing time interpolatation from coarse grid at level "
                             << a_level-1 << endl;
                    }
                    int refratio = get_refratio(domain);

                    // (1) a_curr_fine_step > last_coarse_step : we just advanced the coarse grid
                    // and can now apply time interpolated boundary conditions, or
                    //
                    // (2) a_curr_fine_step < last_coarse_step : we advanced the coarse
                    // grid in a previous step but we still have to do time interpolation (this
                    // can only happen if refratio > 2)

                    exchange_with_coarse_time_interp(domain,a_level,last_coarse_step,
                                                     a_curr_fine_step,refratio);
                    // set_phys_bc(domain,a_level,t_level);
                    a_time_stepper->increment_coarse_exchange_counter(a_level);

                    // Don't increment the fine_exchange_counter, since in the time
                    // interpolated case, there is no coarse time data at time step
                    // a_curr_fine_step a_time_stepper->increment_fine_exchange_counter(a_level-1);
                }
            }
        }
    }
    if (verbose)
    {
        cout << "Taking step on level " << a_level << endl;
    }

    fclaw2d_level_time_data_t time_data;

    time_data.maxcfl = maxcfl;
    time_data.dt = a_time_stepper->dt(a_level);
    time_data.t = t_level;

    // Advance this level from 'a_curr_fine_step' to 'a_curr_fine_step +
    // a_time_stepper.step_inc(a_level)'
    fclaw2d_domain_iterate_level(domain, a_level,
                                 cb_advance_patch,
                                 (void *) &time_data);
    a_time_stepper->increment_step_counter(a_level);
    a_time_stepper->increment_time(a_level);

    level_exchange(domain,a_level);

    set_phys_bc(domain,a_level,t_level);
    a_time_stepper->increment_level_exchange_counter(a_level);

    if (verbose)
    {
        cout << "Advance on level " << a_level << " done" << endl << endl;
    }

    return time_data.maxcfl;  // Maximum from level iteration
}

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
