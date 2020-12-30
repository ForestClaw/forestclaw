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

struct fclaw2d_global;

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
    FCLAW2D_TIMER_NONE = -1,
    FCLAW2D_TIMER_INIT,
    FCLAW2D_TIMER_ADVANCE,
    FCLAW2D_TIMER_ELLIPTIC_SOLVE,
    FCLAW2D_TIMER_GHOSTFILL,
    FCLAW2D_TIMER_REGRID,
    FCLAW2D_TIMER_DIAGNOSTICS,
    FCLAW2D_TIMER_OUTPUT,
    FCLAW2D_TIMER_GHOSTPATCH_COMM,
    FCLAW2D_TIMER_ADAPT_COMM,
    FCLAW2D_TIMER_PARTITION_COMM,
    FCLAW2D_TIMER_DIAGNOSTICS_COMM,
    FCLAW2D_TIMER_CFL_COMM,
    FCLAW2D_TIMER_WALLTIME,
    FCLAW2D_TIMER_UNACCOUNTED,
    FCLAW2D_TIMER_ADVANCE_STEPS_COUNTER,
    FCLAW2D_TIMER_ELLIPTIC_GRIDS_COUNTER,
    FCLAW2D_TIMER_GRIDS_PER_PROC,
    FCLAW2D_TIMER_GRIDS_INTERIOR,
    FCLAW2D_TIMER_GRIDS_LOCAL_BOUNDARY,
    FCLAW2D_TIMER_GRIDS_REMOTE_BOUNDARY,
    FCLAW2D_TIMER_REGRID_BUILD,
    FCLAW2D_TIMER_REGRID_TAGGING,
    FCLAW2D_TIMER_TIMESYNC,
    FCLAW2D_TIMER_GHOSTPATCH_BUILD,
    FCLAW2D_TIMER_PARTITION,
    FCLAW2D_TIMER_PARTITION_BUILD,
    FCLAW2D_TIMER_ADVANCE_STEP2,
    FCLAW2D_TIMER_ADVANCE_B4STEP2,
    FCLAW2D_TIMER_GHOSTFILL_COPY,
    FCLAW2D_TIMER_GHOSTFILL_AVERAGE,
    FCLAW2D_TIMER_GHOSTFILL_INTERP,
    FCLAW2D_TIMER_GHOSTFILL_PHYSBC,
    FCLAW2D_TIMER_GHOSTFILL_STEP1,
    FCLAW2D_TIMER_GHOSTFILL_STEP2,
    FCLAW2D_TIMER_GHOSTFILL_STEP3,
    FCLAW2D_TIMER_NEIGHBOR_SEARCH,
    FCLAW2D_TIMER_LOCAL_COMM,
    FCLAW2D_TIMER_GLOBAL_COMM,
    FCLAW2D_TIMER_CUDA_ALLOCATE,
    FCLAW2D_TIMER_CUDA_MEMCOPY_H2H,
    FCLAW2D_TIMER_CUDA_MEMCOPY_H2D,
    FCLAW2D_TIMER_CUDA_MEMCOPY_D2H,
    FCLAW2D_TIMER_EXTRA1,
    FCLAW2D_TIMER_EXTRA2,
    FCLAW2D_TIMER_EXTRA3,
    FCLAW2D_TIMER_EXTRA4,
    FCLAW2D_TIMER_COUNT
} fclaw2d_timer_names_t;


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
fclaw2d_timer_t;

double fclaw2d_timer_wtime (void);

void fclaw2d_timer_init (fclaw2d_timer_t *timer);

void fclaw2d_timer_start (fclaw2d_timer_t *timer);

void fclaw2d_timer_stop (fclaw2d_timer_t *timer);

void fclaw2d_timer_start_threadsafe(fclaw2d_timer_t *timer);

void fclaw2d_timer_stop_threadsafe(fclaw2d_timer_t *timer);

/* Use keyword 'struct' to avoid circular dependencies */
void fclaw2d_timer_report(struct fclaw2d_global* glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif


#endif
