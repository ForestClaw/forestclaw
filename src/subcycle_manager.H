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

#ifndef SUBCYCLE_MANAGER_H
#define SUBCYCLE_MANAGER_H

#include <iostream>
#include <cstdlib>
#include <vector>

#include "fclaw2d_convenience.h"
#include "forestclaw2d.h"
#include "amr_options.h"

#include "fclaw_defs.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

class level_data
{
public:

    level_data();
    ~level_data();

    void define(const int& a_level,
                const int& a_refratio,
                const int& a_num_patches,
                const int& a_maxlevel,
                const double& a_time,
                const bool& a_subcycle);

    void set_dt(const double& a_dt);

    void increment_step_counter();
    void increment_time();
    double current_time();
    double dt();

    bool level_exchange_done();
    bool exchanged_with_coarser();
    bool exchanged_with_finer();

    void increment_level_exchange_counter();
    void increment_coarse_exchange_counter();
    void increment_fine_exchange_counter();

    int m_level;
    int m_last_step;
    int m_last_level_exchange;
    int m_last_coarse_exchange;
    int m_last_fine_exchange;
    int m_step_inc;
    int m_num_patches;

    double m_time;
    double m_dt;
};


typedef struct fclaw2d_subcycle_info
{
    int coarse_step;
    int fine_step;
    int refratio;
    bool do_time_interp;  // We may not always need to do time interpolation.
    double level_time;
} fclaw2d_subcycle_info_t;


class subcycle_manager
{
public:
    subcycle_manager();
    ~subcycle_manager();
    void define(fclaw2d_domain_t *domain,
                const amr_options_t *gparms,
                const double& a_time);

    bool can_advance(const int& a_level, const int& a_from_step);

    bool solution_updated(const int& a_level, const int& a_step);
    bool level_exchange_done(const int& a_level);
    bool exchanged_with_coarser(const int& a_level);
    bool exchanged_with_finer(const int& a_level);

    int last_step(const int& a_level);
    int step_inc(const int& a_level);
    void increment_step_counter(const int& a_level);

    bool nosubcycle();

    // These deal with real-valued 'time' and 'dt' value.  Most others only deal with integers,
    // i.e. powers of ref_ratio.
    double current_time(const int& a_level);
    double dt(const int& level);
    void increment_time(const int& a_level);
    void set_dt_minlevel(const double& a_dt);

    int minlevel_factor();
    int maxlevel_factor();

    bool is_coarsest(const int& a_level);
    bool is_finest(const int& a_level);
    int minlevel();
    int maxlevel();

    void increment_level_exchange_counter(const int& a_level);
    void increment_coarse_exchange_counter(const int& a_level);
    void increment_fine_exchange_counter(const int& a_level);

    std::vector<level_data> m_levels;
    int m_maxlevel;
    int m_minlevel;
    int m_refratio;
    double m_dt_minlevel;
    double m_t_minlevel;
    bool m_nosubcycle;
    int m_verbosity;
};

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
