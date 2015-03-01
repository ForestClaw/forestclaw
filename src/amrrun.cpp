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


#include <amr_includes.H>
#include <fclaw2d_output.h>
#include <fclaw2d_regrid.h>

#ifdef __cplusplus
extern "C"
{
#if 0
{
#endif
#endif


/*  -----------------------------------------------------------------
    Time stepping
    -- saving time steps
    -- restoring time steps
    -- Time stepping, based on when output files should be created.
    ----------------------------------------------------------------- */

static
void cb_restore_time_step(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    ClawPatch *this_cp = get_clawpatch(this_patch);

    /* Copy most current time step data to grid data. (m_griddata <== m_griddata_last) */
    this_cp->restore_step();
}

static
void restore_time_step(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_iterate_patches(domain,cb_restore_time_step,(void *) NULL);
}

static
void cb_save_time_step(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    ClawPatch *this_cp = get_clawpatch(this_patch);

    /* Copy grid data (m_griddata) on each patch to temporary storage
       (m_griddata_tmp <== m_griddata); */
    this_cp->save_step();
}

static
void save_time_step(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_iterate_patches(domain,cb_save_time_step,(void *) NULL);
}


/* -------------------------------------------------------------------------------
   Output style 1
   Output times are at times [0,dT, 2*dT, 3*dT,...,Tfinal], where dT = tfinal/nout
   -------------------------------------------------------------------------------- */
static void outstyle_1(fclaw2d_domain_t **domain)
{
    int iframe = 0;
    fclaw2d_output_frame(*domain,iframe);

    const amr_options_t *gparms = get_domain_parms(*domain);

    double final_time = gparms->tfinal;
    int nout = gparms->nout;
    double initial_dt = gparms->initial_dt;

    double t0 = 0;

    double dt_outer = (final_time-t0)/double(nout);
    int level_factor = pow_int(2,gparms->minlevel);
    double dt_level0 = initial_dt*level_factor;  // Get a level 0 time step
    double t_curr = t0;

    for(int n = 0; n < nout; n++)
    {
        double tstart = t_curr;
        double tend = tstart + dt_outer;
        int n_inner = 0;
        while (t_curr < tend)
        {
            subcycle_manager time_stepper;
            time_stepper.define(*domain,gparms,t_curr);

            set_domain_time(*domain,t_curr);

            /* In case we have to reject this step */
            if (!gparms->use_fixed_dt)
            {
                save_time_step(*domain);
            }

            if (gparms->run_diagnostics)
            {
                /* Get current domain data since it may change during
                   regrid. */
                fclaw2d_domain_data_t *ddata = get_domain_data(*domain);
                fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_CHECK]);

                fclaw2d_run_diagnostics(*domain);

                fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_CHECK]);
            }

            /* Take a stable level 0 time step (use this as the base
               level time step even if we have no grids on level 0) and
               reduce it. */
            int reduce_factor;
            if (time_stepper.nosubcycle())
            {
                /* Take one step of a stable time step for the finest
                   non-emtpy level. */
                reduce_factor = time_stepper.maxlevel_factor();
            }
            else
            {
                /* Take one step of a stable time step for the coarsest
                   non-empty level. */
                reduce_factor = time_stepper.minlevel_factor();
            }
            double dt_minlevel_desired = dt_level0/reduce_factor;

            /* Use the tolerance to make sure we don't take a tiny time
               step just to hit 'tend'.   We will take a slightly larger
               time step now (dt_cfl + tol) rather than taking a time step
               of 'dt_minlevel' now, followed a time step of only 'tol' in
               the next step.  Of course if 'tend - t_curr > dt_minlevel',
               then dt_minlevel doesn't change. */
            double tol = 1e-2*dt_minlevel_desired;
            fclaw_bool took_small_step = false;
            fclaw_bool took_big_step = false;
            double dt_minlevel = dt_minlevel_desired;
            if (!gparms->use_fixed_dt)
            {
                double small_step = tend-t_curr-dt_minlevel;
                if (small_step  < tol)
                {
                    dt_minlevel = tend - t_curr;  // <= 'dt_minlevel + tol'
                    if (small_step < 0)
                    {
                        /* We have (tend-t_curr) < dt_minlevel, and
                           we have to take a small step to hit tend */
                        took_small_step = true;
                    }
                    else
                    {
                        /* Take a bigger step now to avoid small step
                           in next time step. */
                        took_big_step = true;
                    }
                }
            }
            /* This also sets a scaled time step for all finer levels. */
            time_stepper.set_dt_minlevel(dt_minlevel);

            double maxcfl_step = advance_all_levels(*domain, &time_stepper);

            maxcfl_step = fclaw2d_domain_global_maximum (*domain, maxcfl_step);

            fclaw_global_productionf("Level %d step %5d : dt = %12.3e; maxcfl (step) = " \
                                     "%8.3f; Final time = %12.4f\n",    \
                                     time_stepper.minlevel(),n_inner,dt_minlevel,
                                     maxcfl_step, t_curr+dt_minlevel);

            if (maxcfl_step > gparms->max_cfl)
            {
                fclaw_global_infof("   WARNING : Maximum CFL exceeded; "    \
                                   "retaking time step\n");
                if (!gparms->use_fixed_dt)
                {
                    restore_time_step(*domain);

                    /* Modify dt_level0 from step used. */
                    dt_level0 = dt_level0*gparms->desired_cfl/maxcfl_step;

                    /* Got back to start of loop, without incrementing
                       step counter or time level */
                    continue;
                }
            }

            if (!gparms->use_fixed_dt)
            {
                if (took_small_step)
                {
                    fclaw_global_infof("   WARNING : Took small time step which was " \
                                       "%6.1f%% of desired dt.\n",
                                       100.0*dt_minlevel/dt_minlevel_desired);
                }
                if (took_big_step)
                {
                    fclaw_global_infof("   WARNING : Took big time step which was " \
                                       "%6.1f%% of desired dt.\n",
                                       100.0*dt_minlevel/dt_minlevel_desired);
                }


                /* New time step, which should give a cfl close to the
                   desired cfl. */
                double dt_new = dt_level0*gparms->desired_cfl/maxcfl_step;
                if (!took_small_step)
                {
                    dt_level0 = dt_new;
                }
                else
                {
                    /* use time step that would have been used had we
                       not taken a small step */
                }
            }
            n_inner++;
            t_curr += dt_minlevel;

            if (gparms->regrid_interval > 0)
            {
                if (n_inner % gparms->regrid_interval == 0)
                {
                    fclaw_global_infof("regridding at step %d\n",n);
                    fclaw2d_regrid(domain);
                }
            }
            else
            {
                /* Use a static grid */
            }
        }

        /* Output file at every outer loop iteration */
        set_domain_time(*domain,t_curr);
        iframe++;
        fclaw2d_output_frame(*domain,iframe);
    }
}

#if 0
static void outstyle_2(fclaw2d_domain_t **domain)
{
    // Output time at specific time steps.
}
#endif

static void outstyle_3(fclaw2d_domain_t **domain)
{
    /* Write out an initial time file */
    int iframe = 0;
    fclaw2d_output_frame(*domain,iframe);

    const amr_options_t *gparms = get_domain_parms(*domain);
    double initial_dt = gparms->initial_dt;
    int nstep_outer = gparms->nout;
    int nstep_inner = gparms->nstep;

    double t0 = 0;
    /* The user dt_initial is the appropriate value for minlevel
       (not necessarily level 0) */
    int level_factor = pow_int(2,gparms->minlevel);

    /* Increase dt to value appropriate for level 0 */
    double dt_level0 = initial_dt*level_factor;
    double t_curr = t0;
    set_domain_time(*domain,t_curr);
    int n = 0;

    while (n < nstep_outer)
    {
        int steps_to_minlevel = pow_int(2,(*domain)->global_minlevel - gparms->minlevel);
        for (int m = 0; m < steps_to_minlevel; m++)
        {
            subcycle_manager time_stepper;
            time_stepper.define(*domain,gparms,t_curr);

            /* In case we have to reject this step */
            if (!gparms->use_fixed_dt)
            {
                save_time_step(*domain);
            }

            if (gparms->run_diagnostics)
            {
                /* Get current domain data since it may change during regrid */
                fclaw2d_domain_data_t *ddata = get_domain_data(*domain);
                fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_CHECK]);

                fclaw2d_run_diagnostics(*domain);

                fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_CHECK]);
            }

            /* Take a stable level 0 time step (use level 0 as the
               base level time step even if we have no grids on level 0)
               and reduce it. */
            int reduce_factor;
            if (time_stepper.nosubcycle())
            {
                /* Take one step of a stable time step for the finest
                   non-emtpy level. */
                reduce_factor = time_stepper.maxlevel_factor();
            }
            else
            {
                /* Take one step of a stable time step for the coarsest
                   non-empty level. */
                reduce_factor = time_stepper.minlevel_factor();
            }
            double dt_minlevel = dt_level0/reduce_factor;

            /* This also sets the time step on all finer levels. */
            time_stepper.set_dt_minlevel(dt_minlevel);

            double maxcfl_step = advance_all_levels(*domain, &time_stepper);

            /* This is a collective communication - everybody needs to wait here. */
            maxcfl_step = fclaw2d_domain_global_maximum (*domain, maxcfl_step);

            fclaw_global_productionf("Level %d step %5d : dt = %12.3e; maxcfl (step) = " \
                                     "%8.3f; Final time = %12.4f\n",
                                     time_stepper.minlevel(),n+1,
                                     dt_minlevel,maxcfl_step, t_curr+dt_minlevel);

            if (maxcfl_step > gparms->max_cfl)
            {
                fclaw_global_productionf("   WARNING : Maximum CFL exceeded; retaking time step\n");

                if (!gparms->use_fixed_dt)
                {
                    restore_time_step(*domain);

                    /* Modify dt_level0 from step use */
                    dt_level0 = dt_level0*gparms->desired_cfl/maxcfl_step;

                    /* Got back to start of loop without incrementing step counter or
                       current time. */
                    continue;
                }
            }

            t_curr += dt_minlevel;
            set_domain_time(*domain,t_curr);

            /* New time step, which should give a cfl close to the desired cfl. */
            if (!gparms->use_fixed_dt)
            {
                dt_level0 = dt_level0*gparms->desired_cfl/maxcfl_step;
            }
        }
        n++;

        if (gparms->regrid_interval > 0)
        {
            if (n % gparms->regrid_interval == 0)
            {
                fclaw_global_infof("regridding at step %d\n",n);
                fclaw2d_regrid(domain);
            }
        }
        else
        {
            /* use only initial refinement*/
        }

        if (n % nstep_inner == 0)
        {
            iframe++;
            fclaw2d_output_frame(*domain,iframe);
        }
    }
}

static void outstyle_4(fclaw2d_domain_t **domain)
{
    /* Write out an initial time file */
    int iframe = 0;
    fclaw2d_output_frame(*domain,iframe);

    const amr_options_t *gparms = get_domain_parms(*domain);
    double initial_dt = gparms->initial_dt;
    int nstep_outer = gparms->nout;
    int nstep_inner = gparms->nstep;

    // assume fixed dt that is stable for all grids.

    double t0 = 0;
    double t_curr = t0;
    set_domain_time(*domain,t_curr);
    int n = 0;
    while (n < nstep_outer)
    {
        subcycle_manager time_stepper;
        time_stepper.define(*domain,gparms,t_curr);

        if (gparms->run_diagnostics)
        {
            /* Get current domain data since it may change during regrid */
            fclaw2d_domain_data_t *ddata = get_domain_data(*domain);
            fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_CHECK]);

            fclaw2d_run_diagnostics(*domain);

            fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_CHECK]);
        }

        double dt_minlevel = initial_dt;

        /* This also sets the time step on all finer levels. */
        time_stepper.set_dt_minlevel(dt_minlevel);

        advance_all_levels(*domain, &time_stepper);

        fclaw_global_productionf("Level %d step %5d : dt = %12.3e; Final time = %16.6e\n",
                                 time_stepper.minlevel(),n+1,dt_minlevel, t_curr+dt_minlevel);

        t_curr += dt_minlevel;
        n++;

        set_domain_time(*domain,t_curr);

        if (gparms->regrid_interval > 0)
        {
            if (n % gparms->regrid_interval == 0)
            {
                fclaw_global_infof("regridding at step %d\n",n);

                fclaw2d_regrid(domain);
            }
        }
        else
        {
            /* Only use the initial grid */
        }

        if (n % nstep_inner == 0)
        {
            iframe++;
            fclaw2d_output_frame(*domain,iframe);
        }
    }
}



void amrrun(fclaw2d_domain_t **domain)
{

    const amr_options_t *gparms = get_domain_parms(*domain);

    switch (gparms->outstyle)
    {
    case 1:
        outstyle_1(domain);
        break;
    case 2:
        fclaw_global_essentialf("Outstyle %d not implemented yet\n", gparms->outstyle);
        exit(0);
    case 3:
        outstyle_3(domain);
        break;
    case 4:
        outstyle_4(domain);
        break;
    default:
        fclaw_global_essentialf("Outstyle %d not implemented yet\n", gparms->outstyle);
        exit(0);
    }
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
