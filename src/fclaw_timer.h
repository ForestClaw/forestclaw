/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW_TIMER_H
#define FCLAW_TIMER_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw_global;

/* -----------------------------------------------------------------
   Work with timers
   ----------------------------------------------------------------- */

/* Define struct here to avoid a circular dependency:

    -- fclaw2d_domain.h needs fclaw_timer.h for the definition
    of fclaw2d_timer_t.

    -- fclaw2d_time.h  references an fclaw2d_domain_t struct.
*/

typedef enum
{
    FCLAW_TIMER_NONE = -1,
    FCLAW_TIMER_INIT,
    FCLAW_TIMER_ADVANCE,
    FCLAW_TIMER_ELLIPTIC_SOLVE,
    FCLAW_TIMER_GHOSTFILL,
    FCLAW_TIMER_REGRID,
    FCLAW_TIMER_DIAGNOSTICS,
    FCLAW_TIMER_OUTPUT,
    FCLAW_TIMER_GHOSTPATCH_COMM,
    FCLAW_TIMER_ADAPT_COMM,
    FCLAW_TIMER_PARTITION_COMM,
    FCLAW_TIMER_DIAGNOSTICS_COMM,
    FCLAW_TIMER_CFL_COMM,
    FCLAW_TIMER_WALLTIME,
    FCLAW_TIMER_UNACCOUNTED,
    FCLAW_TIMER_ADVANCE_STEPS_COUNTER,
    FCLAW_TIMER_ELLIPTIC_GRIDS_COUNTER,
    FCLAW_TIMER_GRIDS_PER_PROC,
    FCLAW_TIMER_GRIDS_INTERIOR,
    FCLAW_TIMER_GRIDS_LOCAL_BOUNDARY,
    FCLAW_TIMER_GRIDS_REMOTE_BOUNDARY,
    FCLAW_TIMER_REGRID_BUILD,
    FCLAW_TIMER_REGRID_TAGGING,
    FCLAW_TIMER_TIMESYNC,
    FCLAW_TIMER_GHOSTPATCH_BUILD,
    FCLAW_TIMER_PARTITION,
    FCLAW_TIMER_PARTITION_BUILD,
    FCLAW_TIMER_ADVANCE_STEP2,
    FCLAW_TIMER_ADVANCE_B4STEP2,
    FCLAW_TIMER_GHOSTFILL_COPY,
    FCLAW_TIMER_GHOSTFILL_AVERAGE,
    FCLAW_TIMER_GHOSTFILL_INTERP,
    FCLAW_TIMER_GHOSTFILL_PHYSBC,
    FCLAW_TIMER_GHOSTFILL_STEP1,
    FCLAW_TIMER_GHOSTFILL_STEP2,
    FCLAW_TIMER_GHOSTFILL_STEP3,
    FCLAW_TIMER_NEIGHBOR_SEARCH,
    FCLAW_TIMER_LOCAL_COMM,
    FCLAW_TIMER_GLOBAL_COMM,
    FCLAW_TIMER_CUDA_ALLOCATE,
    FCLAW_TIMER_CUDA_MEMCOPY_H2H,
    FCLAW_TIMER_CUDA_MEMCOPY_H2D,
    FCLAW_TIMER_CUDA_MEMCOPY_D2H,
    FCLAW_TIMER_EXTRA1,
    FCLAW_TIMER_EXTRA2,
    FCLAW_TIMER_EXTRA3,
    FCLAW_TIMER_EXTRA4,
    FCLAW_TIMER_COUNT
} fclaw_timer_names_t;


/* Priority levels for printing timer statistics */
#define  FCLAW_TIMER_PRIORITY_WALL        7
#define  FCLAW_TIMER_PRIORITY_SUMMARY     6
#define  FCLAW_TIMER_PRIORITY_EXCLUSIVE   5
#define  FCLAW_TIMER_PRIORITY_COUNTERS    4
#define  FCLAW_TIMER_PRIORITY_ADVANCE     3
#define  FCLAW_TIMER_PRIORITY_CUDA        2
#define  FCLAW_TIMER_PRIORITY_DETAILS     1
#define  FCLAW_TIMER_PRIORITY_EXTRA       0




typedef struct
{
    int running;
    double started, stopped;
    double cumulative;
}
fclaw_timer_t;

double fclaw_timer_wtime (void);

void fclaw_timer_init (fclaw_timer_t *timer);

void fclaw_timer_start (fclaw_timer_t *timer);

void fclaw_timer_stop (fclaw_timer_t *timer);

void fclaw_timer_start_threadsafe(fclaw_timer_t *timer);

void fclaw_timer_stop_threadsafe(fclaw_timer_t *timer);

/* Use keyword 'struct' to avoid circular dependencies */
void fclaw2d_timer_report(struct fclaw_global* glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif


#endif
