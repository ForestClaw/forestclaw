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

#ifndef FCLAW2D_SUBCYCLE_MANAGER_H
#define FCLAW2D_SUBCYCLE_MANAGER_H

#include <fclaw2d_forestclaw.h>

/* Needed for std::vector<> library */
#include <vector>

class level_data
{
  public:
    level_data() {};
    ~level_data() {};
    int m_last_step;
    int m_step_inc;
    int m_total_steps;  /* Steps this level needs to take */

    double m_time;
    double m_dt;
};

#if 0
typedef struct level_data
{
    int m_last_step;
    int m_step_inc;
    int m_total_steps;  /* Steps this level needs to take */

    double m_time;
    double m_dt;
} level_data_t;
#endif


class subcycle_manager
{
public:
    subcycle_manager();
    ~subcycle_manager();
    void define(fclaw2d_domain_t *domain,
                const amr_options_t *gparms,
                const double time,
                const double dt_step);

    /* Level access functions */
    double level_time(const int level);   /// time() ?
    double dt(int level);
    int step_inc(const int level);
    int steps(int level);

    /* Keep track of time stepping */
    void increment_step_counter(const int level);
    void increment_time(const int level);
    int last_step(const int level);  /* Most recent step taken on this level */

    /* Compute some basic things */
    double compute_alpha(const int level);
    int timeinterp_level(int maxlevel);

private :
    std::vector<level_data> m_levels;
#if 0
    level_data_t *m_levels;
#endif

    int m_local_minlevel;  /* Needed for assertion in computing alpha */

};

#endif
