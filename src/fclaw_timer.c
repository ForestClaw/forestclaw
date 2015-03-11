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

#include <fclaw_timer.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* -----------------------------------------------------------------
   Work with timers
   ----------------------------------------------------------------- */

double
fclaw2d_timer_wtime (void)
{
    return sc_MPI_Wtime ();
}

void
fclaw2d_timer_init (fclaw2d_timer_t *timer)
{
    memset (timer, 0, sizeof (fclaw2d_timer_t));
}

void
fclaw2d_timer_start (fclaw2d_timer_t *timer)
{
    if (!timer->running) {
        timer->started = fclaw2d_timer_wtime ();
        timer->stopped = 0.;
        timer->running = 1;
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void
fclaw2d_timer_stop (fclaw2d_timer_t *timer)
{
    if (timer->running) {
        timer->stopped = fclaw2d_timer_wtime ();
        timer->cumulative += timer->stopped - timer->started;
        timer->running = 0;
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
