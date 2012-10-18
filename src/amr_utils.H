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
#include "amr_options.h"


using namespace std;

int pow_int(int a, int n);

// #define CH_assert(a)

#define CH_SPACEDIM 2
typedef double Real;

const int SpaceDim = CH_SPACEDIM;

const int NumFaces = 4;

class global_parms;

int num_patches(fclaw2d_domain_t *domain, int level);

int pow_int(int a, int n);

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
                     int& maxlevel,
                     int& icycle,
                     int& manifold);

    void setprob_();
}

class level_data
{
public:

    level_data();
    ~level_data();

    void define(const int& a_level,
                const int& a_refratio,
                const int& a_num_patches,
                const int& a_maxlevel,
                const Real& a_time,
                const bool& a_subcycle);

    void set_dt(const Real& a_dt);

    void increment_step_counter();
    void increment_time();
    Real current_time();
    Real dt();

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

    Real m_time;
    Real m_dt;
};

class subcycle_manager
{
public:
    subcycle_manager();
    ~subcycle_manager();
    void define(fclaw2d_domain_t *domain,
                global_parms *gparms,
                const amr_options_t *amropts,
                const Real& a_time);

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
    Real current_time(const int& a_level);
    Real dt(const int& level);
    void increment_time(const int& a_level);
    void set_dt_minlevel(const Real& a_dt);

    int minlevel_factor();
    int maxlevel_factor();

    bool is_coarsest(const int& a_level);
    bool is_finest(const int& a_level);
    int minlevel();
    int maxlevel();

    void increment_level_exchange_counter(const int& a_level);
    void increment_coarse_exchange_counter(const int& a_level);
    void increment_fine_exchange_counter(const int& a_level);

    vector<level_data> m_levels;
    int m_maxlevel;
    int m_minlevel;
    int m_refratio;
    Real m_dt_minlevel;
    Real m_t_minlevel;
    bool m_nosubcycle;
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

    bool m_manifold;
    bool m_mapped;

    // AMR parameters
    int m_refratio;
    int m_minlevel;
    int m_maxlevel;
    bool m_subcycle;

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
    // FArrayBox(const FArrayBox& fbox);
    void define(int a_size, const Box& a_box);
    void define(const Box& a_box, int a_fields);

    Real* dataPtr();
    Box box();
    int size();
    void operator=(const FArrayBox& fbox);
private:
    Real *m_data;
    int m_size;
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
