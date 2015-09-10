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

#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


subcycle_manager::subcycle_manager() {}
subcycle_manager::~subcycle_manager() {}

void subcycle_manager::define(fclaw2d_domain_t *domain,
                              const amr_options_t *gparms,
                              const double a_initial_t)
{
    m_initial_time = a_initial_t; // Time at start of subcycling.
    m_refratio = gparms->refratio;

    /* query the levels that exist in the domain as a whole,
     * regardless of the levels on this processor */
    m_minlevel = domain->global_minlevel;  /* Might be > gparms->minlevel */
    m_maxlevel = domain->global_maxlevel;
    m_levels.resize(m_maxlevel + 1);

    fclaw_bool subcycle = gparms->subcycle;
    m_nosubcycle = !subcycle;
    for (int level = m_minlevel; level <= m_maxlevel; level++)
    {
        m_levels[level].define(level,m_refratio,m_maxlevel,a_initial_t,subcycle);
    }
    m_verbosity = gparms->verbosity;
}

fclaw_bool subcycle_manager::nosubcycle()
{
    return m_nosubcycle;
}

void subcycle_manager::set_dt_minlevel(const double a_dt_minlevel)
{
    // Time step for minimum level (i.e. coarsest non-empty level).
    m_dt_minlevel = a_dt_minlevel;

    double dt_level = a_dt_minlevel;
    m_levels[m_minlevel].set_dt(dt_level);
    for (int level = m_minlevel+1; level <= m_maxlevel; level++)
    {
        if (!nosubcycle())
        {
            dt_level /= m_refratio;
        }
        m_levels[level].set_dt(dt_level);
    }
}


int subcycle_manager::minlevel_factor()
{
    int factor = 1;
    for (int level = 1; level <= m_minlevel; level++)
    {
        factor *= m_refratio;
    }
    return factor;
}

int subcycle_manager::maxlevel_factor()
{
    int factor = 1;
    for (int level = 1; level <= m_maxlevel; level++)
    {
        factor *= m_refratio;
    }
    return factor;
}


int subcycle_manager::minlevel()
{
    return m_minlevel;
}

int subcycle_manager::maxlevel()
{
    return m_maxlevel;
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
    return a_level == m_minlevel;
}

fclaw_bool subcycle_manager::is_finest(const int a_level)
{
    return a_level == m_maxlevel;
}


double subcycle_manager::dt(const int a_level)
{
    return m_levels[a_level].dt();
}

fclaw_bool subcycle_manager::can_advance(const int a_level, const int a_curr_step)
{
    fclaw_bool verbose = false;
    fclaw_bool b1 = solution_updated(a_level,a_curr_step); // do we need this?  We shouldn't be here if
                                                     // we have not taken a time step to 'a_curr_step'
    fclaw_bool b2 = level_exchange_done(a_level);
    fclaw_bool b3 = exchanged_with_coarser(a_level);
    fclaw_bool b4 = exchanged_with_finer(a_level);  // This may not be needed.
    if (verbose)
    {
        if (!b1)
        {
            cout << " --> (can_advance) Solution at level " << a_level <<
                " has not been updated at step " << a_curr_step << endl;
        }
        if (!b2)
        {
            cout << " --> (can_advance) Level exchange at level " << a_level <<
                " has not been done at step " << a_curr_step << endl;
        }
        if (!b3)
        {
            cout << " --> (can_advance) Exchange with coarse grid at level " << a_level-1 <<
                " not done at step " << a_curr_step << endl;
        }
        if (!b4)
        {
            cout << " --> (can_advance) Exchange with finer grid at level " << a_level+1 <<
                " not done at step " << a_curr_step << endl;
        }
    }
    return b1 && b2 && b3 && b4;
}

double subcycle_manager::level_time(const int a_level)
{
    return m_levels[a_level].current_time();
}

double subcycle_manager::initial_time()
{
    return m_initial_time;
}

fclaw_bool subcycle_manager::solution_updated(const int a_level, const int a_step)
{
    return m_levels[a_level].m_last_step >= a_step;
}

fclaw_bool subcycle_manager::level_exchange_done(const int a_level)
{
    return m_levels[a_level].level_exchange_done();
}

fclaw_bool subcycle_manager::exchanged_with_coarser(const int a_level)
{
    if (is_coarsest(a_level))
        return true;
    else
        return m_levels[a_level].exchanged_with_coarser();
}

fclaw_bool subcycle_manager::exchanged_with_finer(const int a_level)
{
    if (is_finest(a_level))
        return true;
    else
        return m_levels[a_level].exchanged_with_finer();
}

void subcycle_manager::increment_level_exchange_counter(const int a_level)
{
    m_levels[a_level].increment_level_exchange_counter();
}

void subcycle_manager::increment_coarse_exchange_counter(const int a_level)
{
    if (!is_coarsest(a_level))
        m_levels[a_level].increment_coarse_exchange_counter();
}

void subcycle_manager::increment_fine_exchange_counter(const int a_level)
{
    if (!is_finest(a_level))
        m_levels[a_level].increment_fine_exchange_counter();
}

int subcycle_manager::step_inc(const int a_level)
{
    return m_levels[a_level].m_step_inc;
}



level_data::level_data() { }
level_data::~level_data() {}

void level_data::define(const int a_level,
                        const int a_refratio,
                        const int a_maxlevel,
                        const double a_time,
                        const fclaw_bool a_subcycle)
{
    m_level = a_level;
    m_last_step = 0;
    m_last_level_exchange = 0;   // Assume that the level exchange has been done at subcycled time
                                 // step 0.
    m_time = a_time;


    // This factor determinines how many finest level steps are equivalent to a single step at
    // this level.
    // Example : 2 level 1 steps are equal to 1 level 0 step, so we have
    //         m_step_inc = 2 for level 0
    //         m_step_inc = 1 for level 1
    m_step_inc = 1;
    if (a_subcycle)
    {
        for (int level = a_maxlevel; level > a_level; level--)
        {
            m_step_inc *= a_refratio;
        }
    }
    m_last_coarse_exchange = -m_step_inc;
    m_last_fine_exchange = -m_step_inc;

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

void level_data::increment_level_exchange_counter()
{
    m_last_level_exchange += m_step_inc;
}

void level_data::increment_coarse_exchange_counter()
{
    m_last_coarse_exchange += m_step_inc;
}

void level_data::increment_fine_exchange_counter()
{
    m_last_fine_exchange += m_step_inc;
}


fclaw_bool level_data::level_exchange_done()
{
    return m_last_level_exchange == m_last_step;
}

fclaw_bool level_data::exchanged_with_coarser()
{
    return m_last_coarse_exchange == m_last_step;
}

fclaw_bool level_data::exchanged_with_finer()
{
    return m_last_fine_exchange == m_last_step;
}
