/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#include <fclaw2d_patch.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw_global.h>
#include <fclaw2d_options.h>
#include <fclaw_advance.h>
#include <fclaw2d_regrid.h>
#include <fclaw2d_output.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_vtable.h>

#include "fclaw_math.h"

/*  -----------------------------------------------------------------
    Time stepping
    -- saving time steps
    -- restoring time steps
    -- Time stepping, based on when output files should be created.
    ----------------------------------------------------------------- */

static
void cb_restore_time_step(fclaw_domain_t *domain,
                          fclaw_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw2d_patch_restore_step(s->glob,this_patch);
}

static
void restore_time_step(fclaw_global_t *glob)
{
    fclaw_global_iterate_patches(glob,cb_restore_time_step,(void *) NULL);

    //fclaw_options_t *fopt = fclaw2d_get_options(glob);
    //fclaw2d_time_sync_reset(glob,fopt->minlevel,fopt->maxlevel,0);
}

static
void cb_save_time_step(fclaw_domain_t *domain,
                       fclaw_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw2d_patch_save_step(s->glob,this_patch);
}

static
void save_time_step(fclaw_global_t *glob)
{
    fclaw_global_iterate_patches(glob,cb_save_time_step,(void *) NULL);
}


/* -------------------------------------------------------------------------------
   Output style 1
   Output times are at times [0,dT, 2*dT, 3*dT,...,Tfinal], where dT = tfinal/nout
   -------------------------------------------------------------------------------- */
static void outstyle_0(fclaw_global_t *glob)
{

    int iframe;

    iframe = 0;
    fclaw2d_output_frame(glob,iframe);

    int init_flag = 1;
    fclaw2d_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    /* Here is where we might include a call to a static solver, for, say,
       an elliptic solve. The initial grid might contain, for example, a
       right hand side. */
#if 0
    vt.time_independent_solve(domain);

    fclaw2d_diagnostics_run(domain);

    iframe++;
    fclaw2d_output_frame(glob,iframe);
#endif

}


/* -------------------------------------------------------------------------------
   Output style 1
   Output times are at times [0,dT, 2*dT, 3*dT,...,Tfinal], where dT = tfinal/nout
   -------------------------------------------------------------------------------- */
static
void outstyle_1(fclaw_global_t *glob)
{
    fclaw_domain_t** domain = &glob->domain;

    /* Set error to 0 */
    int init_flag = 1;  /* Store anything that needs to be stored */
    fclaw2d_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    int iframe = 0;
    fclaw2d_output_frame(glob,iframe);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    double final_time = fclaw_opt->tfinal;
    int nout = fclaw_opt->nout;
    double initial_dt = fclaw_opt->initial_dt;
    int level_factor = pow_int(2,fclaw_opt->maxlevel - fclaw_opt->minlevel);
    double dt_minlevel = initial_dt;

    double t0 = 0;

    double dt_outer = (final_time-t0)/((double) nout);
    double t_curr = t0;
    int n_inner = 0;

    int n;
    for(n = 0; n < nout; n++)
    {
        double tstart = t_curr;

        glob->curr_time = t_curr;
        double tend = tstart + dt_outer;
        while (t_curr < tend)
        {
            /* In case we have to reject this step */
            if (!fclaw_opt->use_fixed_dt)
            {
                save_time_step(glob);
            }

            /* Use the tolerance to make sure we don't take a tiny time
               step just to hit 'tend'.   We will take a slightly larger
               time step now (dt_cfl + tol) rather than taking a time step
               of 'dt_minlevel' now, followed a time step of only 'tol' in
               the next step.  Of course if 'tend - t_curr > dt_minlevel',
               then dt_minlevel doesn't change. */

            double dt_step = dt_minlevel;
            if (fclaw_opt->advance_one_step)
            {
                dt_step /= level_factor;
            }

            double tol = 1e-2*dt_step;
            int took_small_step = 0;
            int took_big_step = 0;
            double dt_step_desired = dt_step;
            if (!fclaw_opt->use_fixed_dt)
            {
                double small_step = tend-(t_curr+dt_step);
                if (small_step  < tol)
                {
                    dt_step = tend - t_curr;  // <= 'dt_minlevel + tol'
                    if (small_step < 0)
                    {
                        /* We have (tend-t_curr) < dt_minlevel, and
                           we have to take a small step to hit tend */
                        took_small_step = 1;
                    }
                    else
                    {
                        /* Take a bigger step now to avoid small step
                           in next time step. */
                        took_big_step = 1;
                    }
                }
            }
            glob->curr_dt = dt_step;  
            double maxcfl_step = fclaw_advance_all_levels(glob, t_curr,dt_step);

            if (fclaw_opt->reduce_cfl)
            {
                /* If we are taking a variable time step, we have to reduce the 
                   maxcfl so that every processor takes the same size dt */
                fclaw_timer_start (&glob->timers[FCLAW_TIMER_CFL_COMM]);
                maxcfl_step = fclaw2d_domain_global_maximum (*domain, maxcfl_step);
                fclaw_timer_stop (&glob->timers[FCLAW_TIMER_CFL_COMM]);                
            }


            double tc = t_curr + dt_step;
            fclaw_global_productionf("Level %d (%d-%d) step %5d : dt = %12.3e; maxcfl (step) = " \
                                     "%16.8f; Final time = %12.4f\n",
                                     fclaw_opt->minlevel,
                                     (*domain)->global_minlevel,
                                     (*domain)->global_maxlevel,
                                     n_inner+1,dt_step,
                                     maxcfl_step, tc);

            if ((maxcfl_step > fclaw_opt->max_cfl) & fclaw_opt->reduce_cfl)
            {
                fclaw_global_essentialf("   WARNING : Maximum CFL exceeded; "    \
                                        "retaking time step\n");

                if (!fclaw_opt->use_fixed_dt)
                {
                    restore_time_step(glob);

                    /* Modify dt_level0 from step used. */
                    dt_minlevel = dt_minlevel*fclaw_opt->desired_cfl/maxcfl_step;

                    /* Got back to start of loop, without incrementing
                       step counter or time level */
                    continue;
                }
            }

            /* We are happy with this step */
            n_inner++;
            t_curr += dt_step;


            /* Update this step, if necessary */
            if (!fclaw_opt->use_fixed_dt)
            {
                double step_fraction = 100.0*dt_step/dt_step_desired;
                if (took_small_step)
                {
                    fclaw_global_infof("   WARNING : Took small time step which was " \
                                       "%6.1f%% of desired dt.\n",
                                       step_fraction);
                }
                if (took_big_step)
                {
                    fclaw_global_infof("   WARNING : Took big time step which was " \
                                       "%6.1f%% of desired dt.\n",step_fraction);

                }


                /* New time step, which should give a cfl close to the
                   desired cfl. */
                double dt_new = dt_minlevel*fclaw_opt->desired_cfl/maxcfl_step;
                if (!took_small_step)
                {
                    dt_minlevel = dt_new;
                }
                else
                {
                    /* use time step that would have been used had we
                       not taken a small step */
                }
            }
            glob->curr_time = t_curr;

            if (fclaw_opt->advance_one_step)
            {
                fclaw2d_diagnostics_gather(glob, init_flag);                
            }

            if (fclaw_opt->regrid_interval > 0)
            {
                if (n_inner % fclaw_opt->regrid_interval == 0)
                {
                    fclaw_global_infof("regridding at step %d\n",n);
                    fclaw2d_regrid(glob);
                }
            }
        }

        /* Output file at every outer loop iteration */
        fclaw2d_diagnostics_gather(glob, init_flag);
        glob->curr_time = t_curr;
        iframe++;
        fclaw2d_output_frame(glob,iframe);
    }
}

#if 0
static void outstyle_2(fclaw_global_t *glob)
{
    // fclaw_domain_t** domain = &glob->domain;
    // Output time at specific time steps.
}
#endif

static
void outstyle_3(fclaw_global_t *glob)
{
    fclaw_domain_t** domain = &glob->domain;

    int init_flag = 1;
    fclaw2d_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    int iframe = 0;
    fclaw2d_output_frame(glob,iframe);


    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double initial_dt = fclaw_opt->initial_dt;


    //fclaw2d_time_sync_reset(glob,fclaw_opt->minlevel,fclaw_opt->maxlevel,1);

    double t0 = 0;
    double dt_minlevel = initial_dt;
    glob->curr_time = t0;
    int nstep_outer = fclaw_opt->nout;
    int nstep_inner = fclaw_opt->nstep;
    int nregrid_interval = fclaw_opt->regrid_interval;
    int level_factor = pow_int(2,fclaw_opt->maxlevel-fclaw_opt->minlevel);
    if (!fclaw_opt->subcycle)
    {
        if (fclaw_opt->advance_one_step)
        {
            if (!fclaw_opt->outstyle_uses_maxlevel)
            {
                /* Multiply nout/nstep by 2^(maxlevel-minlevel) so that
                   a given nout/nstep pair works for both subcycled
                   and non-subcycled cases.
                   Note : Regrid_interval remains unchanged.*/

                nstep_outer *= level_factor;
                nstep_inner *= level_factor;  /* Only produce nout/nstep output files */
            }
        }
    }

    int n = 0;
    double t_curr = t0;
    while (n < nstep_outer)
    {
        double dt_step = dt_minlevel;
        if (!fclaw_opt->subcycle && fclaw_opt->advance_one_step)
        {
            /* if domain->global_maxlevel < fclaw_opt->maxlevel, this choice
               of time step will take more steps than necessary on
               finest level.  */
            dt_step /= level_factor;
        }

        /* In case we have to reject this step */
        if (!fclaw_opt->use_fixed_dt)
        {
            save_time_step(glob);
        }

        /* Get current domain data since it may change during regrid */
        glob->curr_dt = dt_step;
        double maxcfl_step = fclaw_advance_all_levels(glob, t_curr,dt_step);

        /* This is a collective communication - everybody needs to wait here. */
        if (fclaw_opt->reduce_cfl)
        {
            /* If we are taking a variable time step, we have to reduce the 
               maxcfl so that every processor takes the same size dt */
            fclaw_timer_start (&glob->timers[FCLAW_TIMER_CFL_COMM]);
            maxcfl_step = fclaw2d_domain_global_maximum (*domain, maxcfl_step);
            fclaw_timer_stop (&glob->timers[FCLAW_TIMER_CFL_COMM]);     
        }

        double tc = t_curr + dt_step;
        int level2print = (fclaw_opt->advance_one_step && fclaw_opt->outstyle_uses_maxlevel) ?
                          fclaw_opt->maxlevel : fclaw_opt->minlevel;

        fclaw_global_productionf("Level %d (%d-%d) step %5d : dt = %12.3e; maxcfl (step) = " \
                                 "%16.8f; Final time = %12.4f\n",
                                 level2print,
                                 (*domain)->global_minlevel,
                                 (*domain)->global_maxlevel,
                                 n+1,dt_step,maxcfl_step, tc);

        if (fclaw_opt->reduce_cfl & (maxcfl_step > fclaw_opt->max_cfl))
        {
            if (!fclaw_opt->use_fixed_dt)
            {
                fclaw_global_productionf("   WARNING : Maximum CFL exceeded; retaking time step\n");
                restore_time_step(glob);

                dt_minlevel = dt_minlevel*fclaw_opt->desired_cfl/maxcfl_step;

                /* Go back to start of loop without incrementing step counter or
                   current time. */
                continue;
            }
            else
            {
                fclaw_global_productionf("   WARNING : Maximum CFL exceeded\n");
            }
        }

        /* We are happy with this time step */
        t_curr = tc;
        glob->curr_time = t_curr;

        /* New time step, which should give a cfl close to the desired cfl. */
        if (!fclaw_opt->use_fixed_dt)
        {
            dt_minlevel = dt_minlevel*fclaw_opt->desired_cfl/maxcfl_step;
        }

        n++;  /* Increment outer counter */

        if (fclaw_opt->regrid_interval > 0)
        {
            if (n % nregrid_interval == 0)
            {
                fclaw_global_infof("regridding at step %d\n",n);
                fclaw2d_regrid(glob);
            }
        }

        if (fclaw_opt->advance_one_step)
        {
            //fclaw2d_diagnostics_gather(glob,init_flag);
        }

        if (n % nstep_inner == 0)
        {
            iframe++;
            fclaw2d_diagnostics_gather(glob,init_flag);
            fclaw2d_output_frame(glob,iframe);
        }
    }
}


static
void outstyle_4(fclaw_global_t *glob)
{

    /* Write out an initial time file */
    int iframe = 0;
    fclaw2d_output_frame(glob,iframe);

    int init_flag = 1;
    fclaw2d_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double initial_dt = fclaw_opt->initial_dt;
    int nstep_outer = fclaw_opt->nout;
    int nstep_inner = fclaw_opt->nstep;
    double dt_minlevel = initial_dt;

    double t0 = 0;
    double t_curr = t0;
    glob->curr_time = t_curr;
    int n = 0;
    while (n < nstep_outer)
    {
        /* Get current domain data since it may change during regrid */
        fclaw_advance_all_levels(glob, t_curr, dt_minlevel);

        int level2print = (fclaw_opt->advance_one_step && fclaw_opt->outstyle_uses_maxlevel) ?
                          fclaw_opt->maxlevel : fclaw_opt->minlevel;

        fclaw_global_productionf("Level %d step %5d : dt = %12.3e; Final time = %16.6e\n",
                                 level2print,
                                 n+1,dt_minlevel, t_curr+dt_minlevel);

        t_curr += dt_minlevel;
        n++;

        glob->curr_time = t_curr;

        if (fclaw_opt->regrid_interval > 0)
        {
            if (n % fclaw_opt->regrid_interval == 0)
            {
                fclaw_global_infof("regridding at step %d\n",n);

                fclaw2d_regrid(glob);
            }
        }
        else
        {
            /* Only use the initial grid */
        }

        if (n % nstep_inner == 0)
        {
            fclaw2d_diagnostics_gather(glob,init_flag);
            iframe++;
            fclaw2d_output_frame(glob,iframe);
        }
    }
}


/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

void fclaw2d_run(fclaw_global_t *glob)
{

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    switch (fclaw_opt->outstyle)
    {
    case 0:
        outstyle_0(glob);
        break;
    case 1:
        outstyle_1(glob);
        break;
    case 2:
        fclaw_global_essentialf("Outstyle %d not implemented yet\n", fclaw_opt->outstyle);
        exit(0);
    case 3:
        outstyle_3(glob);
        break;
    case 4:
        outstyle_4(glob);
        break;
    default:
        fclaw_global_essentialf("Outstyle %d not implemented yet\n", fclaw_opt->outstyle);
        exit(0);
    }
}
