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
#include <fclaw_math.h>

#include <fclaw2d_subcycle_manager.hpp>

using namespace std;


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


#if 0
subcycle_manager::subcycle_manager()
{
    m_levels = NULL;
}
subcycle_manager::~subcycle_manager()
{
    FCLAW_FREE(m_levels);
}
#endif

subcycle_manager::subcycle_manager() {};
subcycle_manager::~subcycle_manager() {};

void subcycle_manager::define(fclaw2d_domain_t *domain,
                              const amr_options_t *gparms,
                              const double initial_t,
                              const double dt_global_step)
{
    /* Needed for assertion in computing alpha */
    m_local_minlevel = domain->local_minlevel;

    /* Using std::vector<level_data> */
    m_levels.resize(gparms->maxlevel + 1);

#if 0
    m_levels = FCLAW2D_ALLOC(level_data_t,gparms->maxlevel+1);
#endif

    /* Set global indexing for time levels */
    for (int level = gparms->minlevel; level <= gparms->maxlevel; level++)
    {
        m_levels[level].m_last_step = 0;
        m_levels[level].m_time = initial_t;
    }

    /* Set time step and number of steps to take for each level */
    if (!gparms->subcycle)
    {
        int rf = pow_int(2,gparms->maxlevel-gparms->minlevel);
        for (int level = gparms->minlevel; level <= gparms->maxlevel; level++)
        {
            m_levels[level].m_step_inc = 1;
            if (gparms->global_time_stepping)
            {
                /* We do something between each step */
                m_levels[level].m_dt = dt_global_step;
                m_levels[level].m_total_steps = 1;
            }
            else
            {
                /* We take rf steps here without regridding or writing output
                   files between each step.  */
                m_levels[level].m_dt = dt_global_step/rf;
                m_levels[level].m_total_steps = rf;
            }
        }
    }
    else
    {
        double dt_level = dt_global_step;
        int steps_inc = pow_int(2,gparms->maxlevel-gparms->minlevel);
        int total_steps = 1;
        for (int level = gparms->minlevel; level <= gparms->maxlevel; level++)
        {
            m_levels[level].m_dt = dt_level;
            m_levels[level].m_total_steps = total_steps;
            m_levels[level].m_step_inc = steps_inc;
            dt_level /= 2;
            steps_inc /= 2;
            total_steps *= 2;
        }
    }
}

int subcycle_manager::steps(int level)
{
    return m_levels[level].m_total_steps;
}

int subcycle_manager::last_step(const int level)
{
    return m_levels[level].m_last_step;
}

double subcycle_manager::dt(const int level)
{
    return m_levels[level].m_dt;
}

double subcycle_manager::level_time(const int level)
{
    return m_levels[level].m_time;
}


int subcycle_manager::step_inc(const int level)
{
    return m_levels[level].m_step_inc;
}


void subcycle_manager::increment_step_counter(const int level)
{
    int step_inc = m_levels[level].m_step_inc;
    m_levels[level].m_last_step += step_inc;
}

void subcycle_manager::increment_time(const int level)
{
    double dt = m_levels[level].m_dt;
    m_levels[level].m_time += dt;
}

int subcycle_manager::timeinterp_level(int maxlevel)
{
    int ti_level = maxlevel-1;  /* Time interpolated level */
    int last_step = m_levels[maxlevel].m_last_step;
    while (last_step % m_levels[ti_level].m_step_inc == 0)
    {
        ti_level--;
    }
    return ti_level;
}

double subcycle_manager::compute_alpha(const int level)
{
    FCLAW_ASSERT(level > m_local_minlevel);

    /* Time interpolate this data for a future exchange with finer grid */
    int coarse_inc =
        m_levels[level-1].m_step_inc;
    int new_curr_step =
        m_levels[level].m_last_step;
    double alpha =
        double(new_curr_step % coarse_inc)/coarse_inc;

    return alpha;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
