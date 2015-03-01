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

#ifndef FCLAW2D_TYPEDEFS_H
#define FCLAW2D_TYPEDEFS_H

/* this header file must come first */
#include "fclaw2d_defs.H"

#include "fclaw_options.h"
#include "forestclaw2d.h"


class ClawPatch;

typedef struct fclaw2d_level_time_data fclaw2d_level_time_data_t;

typedef void (*fclaw2d_problem_setup_t) (fclaw2d_domain_t * domain);

typedef enum
{
    FCLAW2D_TIMER_NONE = -1,
    FCLAW2D_TIMER_INIT,
    FCLAW2D_TIMER_REGRID,
    FCLAW2D_TIMER_OUTPUT,
    FCLAW2D_TIMER_CHECK,
    FCLAW2D_TIMER_ADVANCE,
    FCLAW2D_TIMER_EXCHANGE,
    FCLAW2D_TIMER_WALLTIME,
    FCLAW2D_TIMER_UNACCOUNTED,
    FCLAW2D_TIMER_BUILDPATCHES,
    FCLAW2D_TIMER_COUNT
}
fclaw2d_timer_names_t;

typedef struct
{
    int running;
    double started, stopped;
    double cumulative;
}
fclaw2d_timer_t;

typedef struct fclaw2d_domain_data
{
    /* Debug counters and timers */
    int count_set_clawpatch, count_delete_clawpatch;
    int count_amr_advance, count_ghost_exchange, count_amr_regrid;
    int is_latest_domain;
    fclaw2d_timer_t timers[FCLAW2D_TIMER_COUNT];

    /* Time at start of each subcycled time step */
    double curr_time;

    /* This should not be copied, but needs to be redone for every new domain */
    fclaw2d_domain_exchange_t *domain_exchange;

}
fclaw2d_domain_data_t;

typedef struct fclaw2d_block_data
{
    int mthbc[FCLAW_NUMFACES];  /* >=0 for physical bc types */
}
fclaw2d_block_data_t;

typedef struct fclaw2d_patch_data
{
    ClawPatch *cp;
}
fclaw2d_patch_data_t;

/* -----------------------------------------------------------
   Data needed for time stepping
   ----------------------------------------------------------- */
struct fclaw2d_level_time_data
{
    /* Single step data. This always has to be set. */
    double dt;
    double t_initial;
    double t_level;
    double t_coarse;

    /* Needed for explicit CFL limited schemes */
    double maxcfl;

    /* Extra data that might be needed for more complicated time stepping.
     * Not always set.
     */
    double alpha;               /* Fraction of coarser dt completed. */
    double dt_coarse;
    bool is_coarsest;
    bool fixed_dt;
};

#endif
