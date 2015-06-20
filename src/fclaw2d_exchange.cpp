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
#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_regrid.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_partition.h>
#include <fclaw2d_exchange.h>

static void
unpack_ghost_patches(fclaw2d_domain_t* domain, fclaw2d_domain_exchange_t *e,
                     int minlevel, int maxlevel, int time_interp)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];
        int level = ghost_patch->level;

        if (level >= minlevel-1)
        {
            int blockno = ghost_patch->u.blockno;

            int patchno = i;

            /* access data stored on remote procs. */
            double *q = (double*) e->ghost_data[patchno];

            int unpack_to_timeinterp_patch=0;
            if (time_interp && level == minlevel-1)
            {
                unpack_to_timeinterp_patch = 1;
            }
            fclaw2d_clawpatch_unpack_ghost(domain, ghost_patch, blockno,
                                           patchno, q, unpack_to_timeinterp_patch);
        }
    }
}

/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* This is called whenever all time levels are time synchronized. */
void fclaw2d_exchange_ghost_patches(fclaw2d_domain_t* domain,
                                    int minlevel,
                                    int maxlevel,
                                    int time_interp)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    fclaw2d_domain_exchange_t *e = fclaw2d_partition_get_exchange_data(domain);

    /* Store pointers to local boundary data.  We do this here
       because we may be exchanging with time interpolated data. */
    int zz = 0;
    for (int nb = 0; nb < domain->num_blocks; ++nb)
    {
        for (int np = 0; np < domain->blocks[nb].num_patches; ++np)
        {
            if (domain->blocks[nb].patches[np].flags &
                FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY)
            {
                fclaw2d_patch_t *this_patch = &domain->blocks[nb].patches[np];
                int level = this_patch->level;
                FCLAW_ASSERT(level <= maxlevel);

                double *q;
                if (time_interp && level == minlevel-1)
                {
                    q = fclaw2d_clawpatch_get_q_timesync(domain,this_patch,time_interp);
                }
                else
                {
                    q = fclaw2d_clawpatch_get_q(domain,this_patch);
                }
                e->patch_data[zz++] = (void*) q;
            }
        }
    }

    /* Do exchange to update ghost patch data */
    if (time_interp)
    {
        int time_interp_level = minlevel-1;
        fclaw2d_domain_ghost_exchange(domain, e, time_interp_level, maxlevel);
        /* Store newly updated e->ghost_patch_data into ghost patches constructed
           locally */
        unpack_ghost_patches(domain,e,time_interp_level,maxlevel,time_interp);
    }
    else
    {
        fclaw2d_domain_ghost_exchange(domain, e, minlevel, maxlevel);
        /* Store newly updated e->ghost_patch_data into ghost patches constructed
           locally */
        unpack_ghost_patches(domain,e,minlevel,maxlevel,time_interp);
    }


    /* Count calls to this function */
    ++ddata->count_ghost_exchange;
}
