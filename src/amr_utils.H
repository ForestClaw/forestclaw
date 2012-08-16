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

#ifndef AMR_UTILS_H
#define AMR_UTILS_H

#include <iostream>
#include <cstdlib>
#include <vector>

#include "fclaw2d_convenience.h"
#include "forestclaw2d.h"

using std::ifstream;
using std::ios;
using std::vector;

#define CH_assert(a)

#define CH_SPACEDIM 2
const int SpaceDim = CH_SPACEDIM;

// Needs to be fixed!
#define Max(a,b)  a

#define cout std::cout

#define endl std::endl

typedef double Real;

class global_parms;

int num_patches(fclaw2d_domain_t *domain, int level);
global_parms* get_domain_data(fclaw2d_domain_t *domain);

extern "C"
{
    void inputparms_(int& mx_leaf,
                     int& my_leaf,
                     Real& initial_dt,
                     Real& tfinal,
                     Real& max_cfl,
                     Real& desired_cfl,
                     int& nout,
                     int& src_term,
                     int& verbose,
                     int& mcapa,
                     int& maux,
                     int& meqn,
                     int& mwaves,
                     int& maxmwaves,
                     int  mthlim[],
                     int& mbc,
                     int  mthbc[],
                     int  order[],
                     int& minlevel,
                     int& maxlevel);

    void setprob_();
}

class level_data
{
public:

    level_data();
    ~level_data();

    void define(const int& a_level,
                const int& a_refratio,
                const int& a_num_patches);

    void set_dt_level(const Real& a_dt);

    int m_level;
    int m_last_time_step;
    int m_last_level_exchange;
    int m_last_coarse_exchange;
    int m_last_fine_exchange;
    Real m_dt_level;
    int m_num_patches;
};

class subcycle_manager
{
public:
    subcycle_manager();
    ~subcycle_manager();
    void define(fclaw2d_domain_t *domain,const Real& a_time);

    bool can_advance(const int& a_level, const int& a_from_step);

    bool time_step_done(const int& a_level, const int& a_time_step);
    int get_last_step(const int& a_level);
    void increment_step_counter(const int& a_level);
    Real current_time(const int& a_level);
    int time_step_inc(const int& a_level);
    void set_dt_minlevel(const Real& a_dt);

    Real reduce_to_minlevel(const Real& a_dt);
    int reduce_to_minlevel_factor();

    bool is_coarsest(const int& a_level);
    bool is_finest(const int& a_level);
    Real get_dt(const int& level);
    int get_minlevel();
    int get_maxlevel();

    bool level_exchange_done(const int& a_level, const int& a_time_step);
    bool coarse_exchange_done(const int& a_level, const int& a_time_step);
    bool fine_exchange_done(const int& a_level, const int& a_time_step);

    void set_level_exchange(const int& a_level, const int& a_time_step);
    void set_coarse_exchange(const int& a_level, const int& a_time_step);
    void set_fine_exchange(const int& a_level, const int& a_time_step);

    vector<level_data> m_levels;
    int m_maxlevel;
    int m_minlevel;
    int m_refratio;
    Real m_dt_minlevel;
    Real m_t_minlevel;
};


class global_parms
{
public:
    int m_mx_leaf;
    int m_my_leaf;

    // Time stepping parameters
    Real m_initial_dt;
    Real m_tfinal;
    Real m_max_cfl;
    Real m_desired_cfl;
    int m_nout;

    // 'method' parameters
#if CH_SPACEDIM == 2
    int m_order[2];
#else
    int m_order[3];
#endif
    int m_verbose;
    int m_src_term;
    int m_mcapa;
    int m_maux;

    int m_method[7];

    // Number of equations, waves and limiters
    int m_meqn;
    int m_mwaves;
    int m_maxmwaves;
    int *m_mthlim;


    int m_mbc;
    int *m_mthbc;

    // AMR parameters
    int m_refratio;
    int m_minlevel;
    int m_maxlevel;

    global_parms();
    ~global_parms();
    void get_inputParams(); // Read from input file
    void print_inputParams();  // Print out for beginning info

};

class Box
{
public:
    Box();
    Box(const int ll[], const int ur[]);
    Box(const Box& a_box);
    int smallEnd(int idir) const;
    int bigEnd(int idir) const;

private:
    int m_ll[2];
    int m_ur[2];
};

class FArrayBox
{
public:

    FArrayBox();
    ~FArrayBox();
    void define(int a_size, const Box& a_box);
    Real* dataPtr() const;
    Box box();
private:
    Real *m_data;
    Box m_box;
};

class ProblemDomain
{
public:
    bool isPeriodic(int idir) const
        {
            return false;
        }
};

#endif
