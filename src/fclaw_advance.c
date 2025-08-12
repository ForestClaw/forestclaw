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

#include <fclaw_advance.h>
#include <fclaw_math.h>

#include <fclaw_timeinterp.h>
#include <fclaw_ghost_fill.h>
#include <fclaw_update_single_step.h>
#include <fclaw_options.h>
#include <fclaw_global.h>
#include <fclaw_vtable.h>

#include <fclaw_time_sync.h>


typedef struct fclaw2d_level_data
{
	int last_step;
	int step_inc;
	int total_steps;  /* This level takes total_steps of step_inc each */

	double initial_time;
	double current_time;
	double dt_step;
} fclaw2d_level_data_t;

/* Rather than over-loading operators ... */
typedef fclaw2d_level_data_t fclaw2d_timestep_counters;

static
void initialize_timestep_counters(fclaw_global_t* glob,
								  fclaw2d_timestep_counters **ts_counter_ptr,
								  double t_init, double dt)
{
	fclaw_domain_t *domain = glob->domain;
	const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
	fclaw2d_timestep_counters *ts_counter;
	int level;

	*ts_counter_ptr = FCLAW_ALLOC(fclaw2d_level_data_t,fclaw_opt->maxlevel+1);

	ts_counter = *ts_counter_ptr;

	/* Set global indexing for time levels */
	for (level = fclaw_opt->minlevel; level <= fclaw_opt->maxlevel; level++)
	{
		ts_counter[level].last_step = 0;
		ts_counter[level].initial_time = t_init;
		ts_counter[level].current_time = t_init;
	}

	/* Set time step and number of steps to take for each level */
	if (fclaw_opt->subcycle)
	{
		double dt_level = dt;
		int steps_inc = pow_int(2,domain->global_maxlevel-fclaw_opt->minlevel);
		int total_steps = 1;
		for (level = fclaw_opt->minlevel; level <= fclaw_opt->maxlevel; level++)
		{
			ts_counter[level].dt_step = dt_level;
			ts_counter[level].total_steps = total_steps;
			ts_counter[level].step_inc = steps_inc;
			dt_level /= 2;
			steps_inc /= 2;
			total_steps *= 2;
		}
	}
	else
	{
		int rf = pow_int(2,fclaw_opt->maxlevel-fclaw_opt->minlevel);
		for (level = fclaw_opt->minlevel; level <= fclaw_opt->maxlevel; level++)
		{
			ts_counter[level].step_inc = 1;
			if (fclaw_opt->advance_one_step)
			{
				/* Take exactly one step of the dt (stable or not!)
				   that was passed into advance;
				   The calling routine wants to do something (regrid,
				   write output, etc) between each step */
				ts_counter[level].dt_step = dt;
				ts_counter[level].total_steps = 1;
			}
			else
			{
				/* Divide the dt passed into into rf smaller steps.  Take
				   this step on all grids.  Return only rf steps have
				   been taken. */
				ts_counter[level].dt_step = dt/rf;
				ts_counter[level].total_steps = rf;
			}
		}
	}
}

static
void delete_timestep_counters(fclaw2d_timestep_counters **ts_counter)
{
	FCLAW_FREE(*ts_counter);
	*ts_counter = NULL;
}

static
void increment_step_counter(fclaw2d_timestep_counters *ts_counter, int level)
{
	int step_inc = ts_counter[level].step_inc;
	ts_counter[level].last_step += step_inc;
}

static
void increment_time(fclaw2d_timestep_counters *ts_counter, int level)
{
	double dt = ts_counter[level].dt_step;
	ts_counter[level].current_time += dt;
}

static
int timeinterp_level(fclaw2d_timestep_counters *ts_counter, int maxlevel)
{
	int ti_level = maxlevel-1;  /* Time interpolated level */
	int last_step = ts_counter[maxlevel].last_step;
	while (last_step % ts_counter[ti_level].step_inc == 0)
	{
		ti_level--;
	}
	return ti_level;
}

static
double compute_alpha(fclaw_global_t *glob,
					 fclaw2d_timestep_counters *ts_counter,
					 int level)
{
	FCLAW_ASSERT(level > glob->domain->local_minlevel);

	/* Time interpolate this data for a future exchange with finer grid */
	int coarse_inc =
		ts_counter[level-1].step_inc;
	int new_curr_step =
		ts_counter[level].last_step;
	double alpha =
		((double) (new_curr_step % coarse_inc))/coarse_inc;

	return alpha;
}



/* ----------------------------------------------------------
   Main time stepping routines
   ---------------------------------------------------------- */

static
double update_level_solution(fclaw_global_t *glob,
							 int level,
							 double t, double dt)
{
	/* There might not be any grids at this level */
	double cfl = fclaw_update_single_step(glob,level,t,dt);

	return cfl;
}

static
double advance_level(fclaw_global_t *glob,
					 const int level,
					 const int curr_fine_step,
					 double maxcfl,
					 fclaw2d_timestep_counters* ts_counter)
{
	fclaw_domain_t* domain = glob->domain;
	const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
	double t_level = ts_counter[level].current_time;
	double dt_level = ts_counter[level].dt_step;

	int this_level = level;
	int coarser_level = level - 1;

	fclaw_global_infof("Advancing level %d from step %d at time %12.6e\n",
					   this_level,curr_fine_step,t_level);
	double cfl_step = update_level_solution(glob,this_level,t_level,dt_level);
	maxcfl = fmax(maxcfl,cfl_step);

	fclaw_global_infof("------ Max CFL on level %d is %12.4e " \
					   " (using dt = %12.4e)\n",this_level,cfl_step,dt_level);

	increment_step_counter(ts_counter,this_level);
	increment_time(ts_counter,this_level);
	if (this_level > domain->local_minlevel)
	{

		int last_coarse_step = ts_counter[coarser_level].last_step;

		if (last_coarse_step == curr_fine_step)
		{
			double cfl_step = advance_level(glob,coarser_level,
											last_coarse_step,
											maxcfl,ts_counter);
			maxcfl = fmax(maxcfl,cfl_step);
			if (fclaw_opt->subcycle)
			{
				double alpha = compute_alpha(glob,ts_counter,this_level);

				fclaw_global_infof("Time interpolating level %d using alpha = %5.2f\n",
								   coarser_level,alpha);

				fclaw_timer_start (&glob->timers[FCLAW_TIMER_EXTRA1]);
				fclaw_timeinterp(glob,coarser_level,alpha);
				fclaw_timer_stop (&glob->timers[FCLAW_TIMER_EXTRA1]);
			}
		}
	}

	fclaw_global_infof("Advance on level %d done at time %12.6e\n\n",
					   level,ts_counter[level].current_time);

	return maxcfl;  /* Maximum from level iteration */
}

/* -------------------------------------------------------------
   Main routine : Called from fclaw2d_run.
   ------------------------------------------------------------- */

double fclaw_advance_all_levels(fclaw_global_t *glob,
								  double t_curr, double dt)
{
	fclaw_domain_t* domain = glob->domain;

	int level;
	fclaw_timer_start (&glob->timers[FCLAW_TIMER_ADVANCE]);

	const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
	fclaw2d_timestep_counters *ts_counter;

	initialize_timestep_counters(glob,&ts_counter,t_curr,dt);

	/* Advance all grids that are present somewhere (on any proc) */
	int minlevel = domain->global_minlevel;
	int maxlevel = domain->global_maxlevel;

	/* Keep track of largest cfl over all grid updates */
	double maxcfl = 0;

	/* Step inc at maxlevel should be 1 by definition */
	FCLAW_ASSERT(ts_counter[maxlevel].step_inc == 1);
	int n_fine_steps = ts_counter[maxlevel].total_steps;
	int nf;
	for(nf = 0; nf < n_fine_steps; nf++)
	{
		/* Coarser levels get updated recursively */
		/* Advance stores anything needed for later synchronization */
		double cfl_step = advance_level(glob,maxlevel,nf,maxcfl,ts_counter);

		maxcfl = fmax(cfl_step,maxcfl);
		int last_step = ts_counter[maxlevel].last_step;

		/* This loop is entered if we are doing local time stepping */
		if (last_step < ts_counter[maxlevel].total_steps)
		{
			/* Do intermediate ghost cell exchanges */
			double sync_time = ts_counter[maxlevel].current_time;
			if (fclaw_opt->subcycle)
			{
				/* Find time interpolated level and do ghost patch exchange
				   and ghost cell exchange for next update. */
				int time_interp_level = timeinterp_level(ts_counter,maxlevel);

				int time_interp = 1;

				/* Do conservative fix up here */
				fclaw_ghost_update(glob,
									 time_interp_level+1,
									 maxlevel,
									 sync_time,
									 time_interp,
									 FCLAW_TIMER_ADVANCE);
				if (fclaw_opt->time_sync)
			    {
			    	fclaw_time_sync(glob,time_interp_level+1,maxlevel);
			    }
			}
			else
			{
				/* End up here is we are doing global time stepping but return from 
				   advance after 2^(maxlevel-minlevel) time steps. */
				int time_interp = 0;
				fclaw_ghost_update(glob,
									 minlevel,
									 maxlevel,
									 sync_time,
									 time_interp,
									 FCLAW_TIMER_ADVANCE);
				if (fclaw_opt->time_sync)
			    {
				    fclaw_time_sync(glob,minlevel,maxlevel);
			    }
			}
		}
	}
	fclaw_global_infof("Advance is done with coarse grid step at " \
					  " time %12.6e\n",t_curr);

	for(level = domain->local_minlevel; level <= domain->global_maxlevel; level++)
	{
		FCLAW_ASSERT(ts_counter[level].last_step ==
					 ts_counter[level].step_inc*ts_counter[level].total_steps);
	}

	double sync_time =  ts_counter[maxlevel].current_time;
	int time_interp = 0;
	fclaw_ghost_update(glob,minlevel,maxlevel,sync_time,
						 time_interp,FCLAW_TIMER_ADVANCE);

	if (fclaw_opt->time_sync)
	{
	    /* This is the final synchronization not covered in fine_steps loop */
		fclaw_time_sync(glob,minlevel,maxlevel);
	}

	delete_timestep_counters(&ts_counter);

	/* Stop the timer */
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_ADVANCE]);

	/* Count total grids on this processor */
	glob->count_grids_per_proc +=  domain->local_num_patches;
	glob->count_grids_local_boundary +=  domain->num_exchange_patches;
	glob->count_grids_remote_boundary +=  domain->num_ghost_patches;

	/* Count the number of times that advance is called */
	++glob->count_amr_advance;

	return maxcfl;
}
