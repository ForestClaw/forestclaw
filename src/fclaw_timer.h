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

#include <p4est_base.h>
#include <fclaw2d_domain.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* -----------------------------------------------------------------
   Work with timers
   ----------------------------------------------------------------- */

#if 0
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
#endif

double fclaw2d_timer_wtime (void);

void fclaw2d_timer_init (fclaw2d_timer_t *timer);

void
    fclaw2d_timer_start (fclaw2d_timer_t *timer);

void
    fclaw2d_timer_stop (fclaw2d_timer_t *timer);

void
    fclaw2d_timer_report(fclaw2d_domain_t* domain);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif


#endif
