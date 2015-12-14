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
                              const double a_initial_t,
                              const double dt_minlevel)
{
    m_initial_time = a_initial_t; // Time at start of subcycling.
    m_refratio = gparms->refratio;

    /* query the levels that exist in the domain as a whole,
     * regardless of the levels on this processor */

    m_local_minlevel = domain->local_minlevel;
    m_local_maxlevel = domain->local_maxlevel;

    m_user_minlevel = gparms->minlevel;
    m_user_maxlevel = gparms->maxlevel;

    m_levels.resize(m_user_maxlevel + 1);

    bool subcycle = gparms->subcycle;
    m_nosubcycle = !subcycle;
    for (int level = m_user_minlevel; level <= m_local_maxlevel; level++)
    {
        m_levels[level].define(level,m_refratio,m_local_maxlevel,a_initial_t,subcycle);
    }
    m_verbosity = gparms->verbosity;
    set_dt_minlevel(dt_minlevel);
}

fclaw_bool subcycle_manager::nosubcycle()
{
    return m_nosubcycle;
}

fclaw_bool subcycle_manager::subcycle()
{
    return !m_nosubcycle;
}

void subcycle_manager::set_dt_minlevel(const double a_dt_minlevel)
{
    double dt_level = a_dt_minlevel;
    int rf = pow(m_refratio,m_user_maxlevel-m_user_minlevel);
    if (!subcycle())
    {
        dt_level /= rf;
    }
    m_dt_minlevel = dt_level;

    m_levels[m_user_minlevel].set_dt(dt_level);
    for (int level = m_user_minlevel+1; level <= m_user_maxlevel; level++)
    {
        if (subcycle())
        {
            dt_level /= m_refratio;
        }
        m_levels[level].set_dt(dt_level);
    }
}

void subcycle_manager::set_dt_maxlevel(const double a_dt_maxlevel)
{
    double dt_level = a_dt_maxlevel;
    m_levels[m_user_maxlevel].set_dt(dt_level);
    for (int level = m_user_minlevel; level <= m_user_maxlevel; level++)
    {
        m_levels[level].set_dt(dt_level);
    }
}


int subcycle_manager::minlevel_factor()
{
    int factor = 1;
    for (int level = 1; level <= m_user_minlevel; level++)
    {
        factor *= m_refratio;
    }
    return factor;
}

int subcycle_manager::maxlevel_factor()
{
    int factor = 1;
    for (int level = 1; level <= m_user_maxlevel; level++)
    {
        factor *= m_refratio;
    }
    return factor;
}

int subcycle_manager::minlevel()
{
    return m_user_minlevel;
}

int subcycle_manager::maxlevel()
{
    return m_user_maxlevel;
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



int subcycle_manager::verbosity()
{
    return m_verbosity;
}

int subcycle_manager::last_step(const int a_level)
{
    return m_levels[a_level].m_last_step;
}

void subcycle_manager::increment_step_counter(const int a_level)
{
    m_levels[a_level].increment_step_counter();
}

void subcycle_manager::increment_time(const int a_level)
{
    m_levels[a_level].increment_time();
}

fclaw_bool subcycle_manager::is_coarsest(const int a_level)
{
    return a_level == m_local_minlevel;
}


double subcycle_manager::dt(const int a_level)
{
return m_levels[a_level].dt();
}


double subcycle_manager::level_time(const int a_level)
{
    return m_levels[a_level].current_time();
}

double subcycle_manager::sync_time()
{
    return m_levels[m_local_maxlevel].current_time();
}

double subcycle_manager::initial_time()
{
    return m_initial_time;
}


int subcycle_manager::step_inc(const int a_level)
{
    int step_inc = m_levels[a_level].m_step_inc;
    if (subcycle())
    {
        /* step_inc /= m_levels[m_local_maxlevel].m_step_inc; */
    }
    return step_inc;
}



level_data::level_data() { }
level_data::~level_data() {}

void level_data::define(const int level,
                        const int refratio,
                        const int maxlevel,
                        const double time,
                        const fclaw_bool subcycle)
{
    m_level = level;
    m_last_step = 0;
    m_time = time;

    if (subcycle)
    {
        m_step_inc = pow_int(refratio,maxlevel-level);
    }
    else
    {
        m_step_inc = 1;
    }
}

void level_data::set_dt(const double a_dt_level)
{
    m_dt = a_dt_level;
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

void level_data::set_time(const double a_t_level)
{
    m_time = a_t_level;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
