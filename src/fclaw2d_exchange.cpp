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

static void
unpack_ghost_patches(fclaw2d_domain_t* domain, fclaw2d_domain_exchange_t *e,
                     int exchange_minlevel, int exchange_maxlevel)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];
        if (exchange_minlevel <= ghost_patch->level &&
            ghost_patch->level <= exchange_maxlevel)
        {
            int blockno = ghost_patch->u.blockno;

            int patchno = i;

            /* access data stored on remote procs. */
            double *q = (double*) e->ghost_data[patchno];

            fclaw_bool time_interp = ghost_patch->level == exchange_minlevel;
            fclaw2d_clawpatch_unpack_ghost(domain, ghost_patch,blockno,
                                           patchno, q, time_interp);
        }
    }
}

static void
unpack_ghost_patches_all(fclaw2d_domain_t* domain, fclaw2d_domain_exchange_t *e)
{
    // fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];
        int blockno = ghost_patch->u.blockno;

        int patchno = i;

        /* access data stored on remote procs.  */
        double *q = (double*) e->ghost_data[patchno];

        fclaw_bool time_interp = fclaw_false;
        fclaw2d_clawpatch_unpack_ghost(domain, ghost_patch,blockno,
                                       patchno, q, time_interp);
    }
}

/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */

/* This is called whenever all time levels are time synchronized. */
void fclaw2d_exchange_ghost_patches_all(fclaw2d_domain_t* domain)
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
                double *q = fclaw2d_clawpatch_get_q(domain,this_patch);
                e->patch_data[zz++] = (void*) q;
            }
        }
    }

    int minlevel = domain->global_minlevel;
    int maxlevel = domain->global_maxlevel;

    /* Do exchange to update ghost patch data */
    fclaw2d_domain_ghost_exchange(domain, e, minlevel, maxlevel);

    /* Store newly updated e->ghost_patch_data into ghost patches constructed
       locally */
    unpack_ghost_patches_all(domain,e);

    /* Count calls to this function */
    ++ddata->count_ghost_exchange;
}


/* ------------------------------------------------------------------
   Partial exchange.  This differs from the "exchange all" routines
   above because here, we can always work with existing ghost patch
   buffers;  we are not building new ones.
 -------------------------------------------------------------------- */
/* Exchange_minlevel is a time interpolated level. */
void fclaw2d_exchange_stage_data(fclaw2d_domain_t* domain,
                                 int exchange_minlevel,
                                 int exchange_maxlevel)
{
    // fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    fclaw2d_domain_exchange_t *e = fclaw2d_partition_get_exchange_data(domain);

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
                FCLAW_ASSERT(level <= exchange_maxlevel);

                double *q;
                int time_interp = level == exchange_minlevel;
                q = fclaw2d_clawpatch_get_q_timesync(domain,this_patch,time_interp);
                e->patch_data[zz++] = (void*) q;        /* Put this patch's data location */
            }
        }
    }
}

/* This is called anytime we need to update ghost patch data for certain levels
   The assumption is that the finest level is a time_interpolated level.  The
   routine 'unpack_ghost_patches' knows this, and so unpacks ghost patches to the
   correct places.
 */
void fclaw2d_exchange_ghost_patch_partial(fclaw2d_domain_t* domain,
                                          int exchange_minlevel,
                                          int exchange_maxlevel)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    fclaw2d_domain_exchange_t *e = fclaw2d_partition_get_exchange_data(domain);

    /* Do exchange to update ghost patch data */
    fclaw2d_domain_ghost_exchange(domain, e,
                                  exchange_minlevel, exchange_maxlevel);

    /* Store newly updated e->ghost_patch_data into ghost patches constructed
       locally */
    unpack_ghost_patches(domain,e, exchange_minlevel, exchange_maxlevel);

    /* Count calls to this function */
    ++ddata->count_ghost_exchange;
}
