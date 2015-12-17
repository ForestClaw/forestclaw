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
#include <fclaw2d_subcycle_manager.hpp>
#include <fclaw_math.h>

#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif



subcycle_manager::subcycle_manager() {}
subcycle_manager::~subcycle_manager() {}

void subcycle_manager::define(fclaw2d_domain_t *domain,
                              const amr_options_t *gparms,
                              const double initial_t,
                              const double dt_global_step)
{
    m_refratio = gparms->refratio;
    m_initial_time = initial_t;
    m_local_minlevel = domain->local_minlevel;
    m_local_maxlevel = domain->local_maxlevel;

    m_global_minlevel = domain->global_minlevel;
    m_global_maxlevel = domain->global_maxlevel;

    m_user_minlevel = gparms->minlevel;
    m_user_maxlevel = gparms->maxlevel;

    m_levels.resize(m_user_maxlevel + 1);

    m_subcycle = gparms->subcycle;
    m_global_time_stepping = gparms->global_time_stepping;

    /* Set global indexing for time levels */
    for (int level = m_user_minlevel; level <= m_user_maxlevel; level++)
    {
        m_levels[level].define(level,initial_t);
    }

    /* Set time step and number of steps to take for each level */
    m_global_step = dt_global_step;
    if (!m_subcycle)
    {
        int rf = pow_int(2,m_user_maxlevel-m_user_minlevel);
        for (int level = m_user_minlevel; level <= m_user_maxlevel; level++)
        {
            m_levels[level].m_step_inc = 1;
            if (m_global_time_stepping)
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
        int steps_inc = pow_int(2,m_user_maxlevel-m_user_minlevel);
        int total_steps = 1;
        for (int level = m_user_minlevel; level <= m_user_maxlevel; level++)
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

int subcycle_manager::global_time_stepping()
{
    return m_global_time_stepping;
}

double subcycle_manager::global_step()
{
    return m_global_step;
}

fclaw_bool subcycle_manager::subcycle()
{
    return m_subcycle;
}


int subcycle_manager::fine_steps()
{
    return m_levels[m_global_maxlevel].m_total_steps;
}

int subcycle_manager::local_minlevel()
{
    return m_local_minlevel;
}

int subcycle_manager::local_maxlevel()
{
    return m_local_maxlevel;
}

int subcycle_manager::user_minlevel()
{
    return m_user_minlevel;
}


int subcycle_manager::user_maxlevel()
{
    return m_user_maxlevel;
}

int subcycle_manager::global_minlevel()
{
    return m_global_minlevel;
}


int subcycle_manager::global_maxlevel()
{
    return m_global_maxlevel;
}



int subcycle_manager::last_step(const int a_level)
{
    return m_levels[a_level].m_last_step;
}

#if 0
int subcycle_manager::is_final_step(const int level)
{
    return m_levels[level].m_last_step == m_levels[level].m_total_steps;
}
#endif


void subcycle_manager::increment_step_counter(const int level)
{
    m_levels[level].increment_step_counter();
}

void subcycle_manager::increment_time(const int level)
{
    m_levels[level].increment_time();
}

fclaw_bool subcycle_manager::is_coarsest(const int level)
{
    return level == m_local_minlevel;
}


double subcycle_manager::dt(const int level)
{
    return m_levels[level].dt();
}


double subcycle_manager::level_time(const int level)
{
    return m_levels[level].current_time();
}

double subcycle_manager::sync_time()
{
    return m_levels[m_local_maxlevel].current_time();
}

double subcycle_manager::initial_time()
{
    return m_initial_time;
}


int subcycle_manager::step_inc(const int level)
{
    return m_levels[level].m_step_inc;
}



level_data::level_data() {}
level_data::~level_data() {}

void level_data::define(const int level,
                        const double time)
{
    m_level = level;
    m_last_step = 0;
    m_time = time;
}

void level_data::set_dt(const double dt_level)
{
    m_dt = dt_level;
}

double level_data::dt()
{
    return m_dt;
}


void level_data::increment_step_counter()
{
    m_last_step += m_step_inc;
}

void level_data::increment_time()
{
    m_time += m_dt;
}

double level_data::current_time()
{
    return m_time;
}

void level_data::set_time(const double t_level)
{
    m_time = t_level;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
