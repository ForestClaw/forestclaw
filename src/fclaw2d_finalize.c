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

#include <sc_statistics.h>

/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

void fclaw2d_finalize(fclaw2d_domain_t **domain)
{
    const amr_options_t *gparms = get_domain_parms(*domain);

    fclaw_global_essentialf("Finalizing run\n");
    fclaw2d_domain_barrier (*domain);

    if (gparms->report_timing)
    {
        if (gparms->outstyle > 0)
        {
            /* Only call this if we have taken time steps.  For time-independent problems, we
               probably need a different report (no "amr_advance_steps") */
            fclaw2d_timer_report(*domain);
        }
        else
        {
            fclaw_global_essentialf("Timing reports not generated for outstyle=0\n");
        }
    }

    fclaw2d_domain_reset(domain);
}
