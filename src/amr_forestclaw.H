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

#ifndef AMR_FORESTCLAW_H
#define AMR_FORESTCLAW_H

#include "forestclaw2d.h"
#include "amr_options.h"
#include "amr_utils.H"
#include "ClawPatch.H"

using namespace std;

class ClawPatch;

// using std::vector;

typedef struct fclaw2d_domain_data
{
    global_parms *parms;
    const amr_options_t *amropts;
    double curr_time;
} fclaw2d_domain_data_t;

typedef struct fclaw2d_block_data
{
    int mthbc[4];               /* >0 for physical bc types */
} fclaw2d_block_data_t;

typedef struct fclaw2d_patch_data
{
    ClawPatch	*cp;
} fclaw2d_patch_data_t;

typedef struct fclaw2d_level_time_data
{
    Real dt;
    Real maxcfl;
    Real t;
} fclaw2d_level_time_data_t;

typedef struct fclaw2d_subcycle_info
{
    int coarse_step;
    int fine_step;
    int refratio;
    bool do_time_interp;  // We may not always need to do time interpolation.
} fclaw2d_subcycle_info_t;


extern "C"
{
    void write_tfile_(int *iframe,
                      double *time,
                      int *meqn,
                      int *ngrids,
                      int *maux);

    void new_qfile_(int* iframe);

}

fclaw2d_domain_data_t *get_domain_data(fclaw2d_domain_t *domain);
void set_domain_data(fclaw2d_domain_t *domain, global_parms *parms);
void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp);
ClawPatch* get_patch_data(fclaw2d_patch_t *patch);

double get_domain_time(fclaw2d_domain_t *domain);
void set_domain_time(fclaw2d_domain_t *domain,double t);

fclaw2d_block_data_t *get_block_data(fclaw2d_block_t *block);

void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp);

void bc_level_exchange(fclaw2d_domain_t *domain, const int& a_level);

void initial_grid(fclaw2d_domain_t *domain);

void amrsetup(fclaw2d_domain_t *domain);
void amrinit(fclaw2d_domain_t *domain,
             global_parms *parms, const amr_options_t * amropts);
void amrrun(fclaw2d_domain_t *domain);
void amrout(fclaw2d_domain_t *domain, int iframe);
void amrreset(fclaw2d_domain_t *domain);

#endif
