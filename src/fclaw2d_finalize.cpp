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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_partition.h>

#include <sc_statistics.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


#if 0
#define FCLAW2D_STATS_SET(stats,ddata,NAME) do {                               \
    SC_CHECK_ABORT (!(ddata)->timers[FCLAW2D_TIMER_ ## NAME].running,          \
                    "Timer " #NAME " still running in amrreset");              \
    sc_stats_set1 ((stats) + FCLAW2D_TIMER_ ## NAME,                           \
                   (ddata)->timers[FCLAW2D_TIMER_ ## NAME].cumulative, #NAME); \
} while (0)
#endif


/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

void fclaw2d_finalize(fclaw2d_domain_t **domain)
{

    fclaw_global_essentialf("Finalizing run\n");
    fclaw2d_domain_barrier (*domain);

    fclaw2d_timer_report(*domain);

    fclaw2d_domain_reset(domain);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
