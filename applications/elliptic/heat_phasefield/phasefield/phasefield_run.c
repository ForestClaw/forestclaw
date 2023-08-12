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

#include "phasefield_user.h"
#include "phasefield_operator.h"

#include <fclaw_clawpatch.h>
#include <fclaw_patch.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw_global.h>
#include <fclaw_ghost_fill.h>
#include <fclaw_options.h>
#include <fclaw_advance.h>
#include <fclaw_regrid.h>
#include <fclaw_output.h>
#include <fclaw_diagnostics.h>
#include <fclaw_vtable.h>

#include <fclaw_elliptic_solver.h>

#include "fclaw_math.h"


/*  -----------------------------------------------------------------
    Time stepping
    -- saving time steps
    -- restoring time steps
    -- Time stepping, based on when output files should be created.
    ----------------------------------------------------------------- */

static
void update_q(fclaw_domain_t *domain,
              fclaw_patch_t *patch,
              int blockno, int patchno,
              void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t *) user;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(g->glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double* rhs;
    int mfields;
    fclaw_clawpatch_rhs_data(g->glob, patch, &rhs, &mfields);

    double* q;
    int meqn;
    fclaw_clawpatch_soln_data(g->glob, patch, &q, &meqn);

    FCLAW_ASSERT(mfields==meqn);

    PHASEFIELD_UPDATE_Q(&mbc,&mx,&my,&meqn,&mfields,rhs,q);
} 

static
void phasefield_run_update_q(fclaw_global_t *glob)
{
    fclaw_global_iterate_patches(glob, update_q, NULL);
}


/* -------------------------------------------------------------------------------
   Output style 1
   Output times are at times [0,dT, 2*dT, 3*dT,...,Tfinal], where dT = tfinal/nout
   -------------------------------------------------------------------------------- */
static
void outstyle_1(fclaw_global_t *glob)
{
    fclaw_domain_t** domain = &glob->domain;

    int iframe = 0;

    fclaw_output_frame(glob,iframe);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    double final_time = fclaw_opt->tfinal;
    int nout = fclaw_opt->nout;
    double initial_dt = fclaw_opt->initial_dt;
    int level_factor = pow_int(2,fclaw_opt->maxlevel - fclaw_opt->minlevel);
    double dt_minlevel = initial_dt;


    int init_flag = 1;  /* Store anything that needs to be stored */
    fclaw_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    double t0 = 0;

    double dt_outer = (final_time-t0)/((double) nout);
    double t_curr = t0;
    int n_inner = 0;
    //double dt_step_fixed = fclaw_opt->initial_dt;


    int n;
    for(n = 0; n < nout; n++)
    {
        double tstart = t_curr;

        glob->curr_time = t_curr;
        double tend = tstart + dt_outer;
        while (t_curr < tend)
        {
            double dt_step = dt_minlevel;
            if (fclaw_opt->advance_one_step)
            {
                dt_step /= level_factor;
            }

            double tol = 0.01;
            double dt_small = t_curr + dt_step - tend;
            if (dt_small < 0)
            {
                /* May end up taking either larger or smaller step */
                if (fabs(dt_small) < tol*dt_step)
                    dt_step = tend - t_curr;
                //double p = dt_step/dt_step_fixed;
                //fclaw_global_essentialf("Took modified time step; dt = %.4f (%.2f\%)\n",dt_step,p);
            }
            else
            {
                dt_step = tend - t_curr;
            }

            glob->curr_dt = dt_step;  

            /* Set lambda for Backward Euler */
            double lambda = -1/dt_step;
            phasefield_set_lambda(lambda);

            /* Solve the elliptic problem; RHS is set here */
            fclaw_elliptic_solve(glob);

            /* Update solution stored in RHS */
            phasefield_run_update_q(glob);

            double tc = t_curr + dt_step;
            fclaw_global_productionf("Level %d (%d-%d) step %5d : dt = %12.3e; " \
                                     "Final time = %12.4f\n",
                                     fclaw_opt->minlevel,
                                     (*domain)->global_minlevel,
                                     (*domain)->global_maxlevel,
                                     n_inner+1,dt_step,
                                     tc);

            n_inner++;
            t_curr += dt_step;

            glob->curr_time = t_curr;

            fclaw_diagnostics_gather(glob, init_flag);                

            if (fclaw_opt->regrid_interval > 0)
                if (n_inner % fclaw_opt->regrid_interval == 0)
                {
                    fclaw_global_infof("regridding at step %d\n",n);
                    fclaw_regrid(glob);
                }
        }

        /* Output file at every outer loop iteration */
        fclaw_diagnostics_gather(glob, init_flag);
        glob->curr_time = t_curr;
        iframe++;
        fclaw_output_frame(glob,iframe);
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
    fclaw_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    int iframe = 0;
    fclaw_output_frame(glob,iframe);


    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
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

        /* Get current domain data since it may change during regrid */
        glob->curr_dt = dt_step;
        /* Set lambda for Backward Euler */
        double lambda = -1/dt_step;
        phasefield_set_lambda(lambda);


        /* Solve the elliptic problem; RHS is set here */
        fclaw_elliptic_solve(glob);

        /* Update solution stored in RHS */
        phasefield_run_update_q(glob);

        int time_interp = 0;
        fclaw_ghost_update(glob,fclaw_opt->minlevel,fclaw_opt->maxlevel,t_curr,
                             time_interp,FCLAW_TIMER_NONE);


        double tc = t_curr + dt_step;
        int level2print = (fclaw_opt->advance_one_step && fclaw_opt->outstyle_uses_maxlevel) ?
                          fclaw_opt->maxlevel : fclaw_opt->minlevel;

        fclaw_global_productionf("Level %d (%d-%d) step %5d : dt = %12.3e; " \
                                 "Final time = %12.4f\n",
                                 level2print,
                                 (*domain)->global_minlevel,
                                 (*domain)->global_maxlevel,
                                 n+1,dt_step,tc);

        t_curr = tc;
        glob->curr_time = t_curr;

        n++;  /* Increment outer counter */

        if (fclaw_opt->regrid_interval > 0)
        {
            if (n % nregrid_interval == 0)
            {
                fclaw_global_infof("regridding at step %d\n",n);
                fclaw_regrid(glob);
            }
        }

        if (fclaw_opt->advance_one_step)
            fclaw_diagnostics_gather(glob,init_flag);
        

        if (n % nstep_inner == 0)
        {
            iframe++;
            //fclaw2d_diagnostics_gather(glob,init_flag);
            fclaw_output_frame(glob,iframe);
        }
    }
}


static
void outstyle_4(fclaw_global_t *glob)
{

    /* Write out an initial time file */
    int iframe = 0;
    fclaw_output_frame(glob,iframe);

    int init_flag = 1;
    fclaw_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
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

                fclaw_regrid(glob);
            }
        }
        else
        {
            /* Only use the initial grid */
        }

        if (n % nstep_inner == 0)
        {
            fclaw_diagnostics_gather(glob,init_flag);
            iframe++;
            fclaw_output_frame(glob,iframe);
        }
    }
}


/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

void phasefield_run(fclaw_global_t *glob)
{

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    switch (fclaw_opt->outstyle)
    {
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
