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

#include <fclaw2d_timeinterp.h>
#include <fclaw_patch.h>
#include <fclaw_global.h>

static
void cb_setup_time_interp(fclaw_domain_t *domain,
                          fclaw_patch_t *this_patch,
                          int blockno,
                          int patchno,
                          void *user)
{
    fclaw_global_iterate_t *s = (fclaw_global_iterate_t*) user;
    if (fclaw_patch_has_finegrid_neighbors(this_patch))
    {
        double alpha = *((double*) s->user);
        fclaw_patch_setup_timeinterp(s->glob,this_patch,alpha);
    }
}

/* ----------------------------------------------------------------------
   Main routine in this file.  This file assumes that both coarse and
   fine grids have valid interior data;  they only need to exchange (
   via interpolating and averaging) ghost cell values.
   -------------------------------------------------------------------- */

void fclaw2d_timeinterp(fclaw_global_t *glob,
                        int level,double alpha)
{
    /* Store time interpolated data into m_griddata_time_sync. */
    fclaw_global_iterate_level(glob,level,cb_setup_time_interp,
                                 (void *) &alpha);
}
