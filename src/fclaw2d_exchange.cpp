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
#include <fclaw2d_domain.h>
#include <fclaw2d_regrid.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_partition.h>
#include <fclaw2d_exchange.h>


/* Also needed in amrreset */
fclaw2d_domain_exchange_t*
    fclaw2d_exchange_get_data(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (domain);
    return ddata->domain_exchange;
}

static
void set_exchange_data(fclaw2d_domain_t* domain,
                       fclaw2d_domain_exchange_t *e)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (domain);
    ddata->domain_exchange = e;
}

static
void build_ghost_patches(fclaw2d_domain_t* domain)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];

        int blockno = ghost_patch->u.blockno;

        /* not clear how useful this patchno is.  In any case, it isn't
           used in defining the ClawPatch, so probably doesn't
           need to be passed in */
        int patchno = i;

        fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_GHOST;

        fclaw2d_patch_data_new(domain,ghost_patch);
        fclaw2d_clawpatch_build_cb(domain,ghost_patch,blockno,
                                   patchno,(void*) &build_mode);
    }
}

static
void delete_ghost_patches(fclaw2d_domain_t* domain)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];
        fclaw2d_patch_data_delete(domain,ghost_patch);
    }
}


static void
unpack_ghost_patches(fclaw2d_domain_t* domain,
                     fclaw2d_domain_exchange_t *e,
                     int minlevel,
                     int maxlevel,
                     int time_interp)
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
            fclaw2d_clawpatch_ghost_unpack(domain, ghost_patch, blockno,
                                           patchno, q, unpack_to_timeinterp_patch);
        }
    }
}

/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* This is called by rebuild_domain */
void fclaw2d_exchange_setup(fclaw2d_domain* domain)
{
    size_t data_size =  fclaw2d_clawpatch_ghost_packsize(domain);
    fclaw2d_domain_exchange_t *e;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(domain);


    /* we just created a grid by amrinit or regrid and we now need to
       allocate data to store and retrieve local boundary patches and
       remote ghost patches */
    e = fclaw2d_domain_allocate_before_exchange (domain, data_size);

    int zz = 0;
    for (int nb = 0; nb < domain->num_blocks; ++nb)
    {
        for (int np = 0; np < domain->blocks[nb].num_patches; ++np)
        {
            if (domain->blocks[nb].patches[np].flags &
                FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY)
            {
                /* Copy q and area into one contingous block */
                fclaw2d_clawpatch_ghost_pack_location(domain,this_patch,
                                                      &e->patch_data[zz++]);
            }
        }
    }


    /* Store e so we can retrieve it later */
    set_exchange_data(domain,e);

    /* Build patches that can be filled later with q data */
    build_ghost_patches(domain);
}

void fclaw2d_exchange_delete(fclaw2d_domain_t** domain)
{
    /* Free old parallel ghost patch data structure, must exist by construction. */
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (*domain);
    fclaw2d_domain_exchange_t *e_old = ddata->domain_exchange;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(*domain);

    if (gparms->manifold)
    {
        if (e_old != NULL)
        {
            /* Delete local patches which are passed to other procs */
            int zz = 0;
            for (int nb = 0; nb < (*domain)->num_blocks; ++nb)
            {
                for (int np = 0; np < (*domain)->blocks[nb].num_patches; ++np)
                {
                    if ((*domain)->blocks[nb].patches[np].flags &
                        FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY)
                    {
                        FCLAW_FREE(e_old->patch_data[zz++]);
                    }
                }
            }
        }
    }
    /* Delete ghost patches */
    delete_ghost_patches(*domain);
    fclaw2d_domain_free_after_exchange (*domain, e_old);
}


/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */

/* This is called whenever all time levels are time synchronized. */
void fclaw2d_exchange_ghost_patches(fclaw2d_domain_t* domain,
                                    int minlevel,
                                    int maxlevel,
                                    int time_interp)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (domain);
    fclaw2d_domain_exchange_t *e = fclaw2d_exchange_get_data(domain);
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(domain);


    /* Store pointers to local boundary data.  We do this here
       because we may be exchanging with time interpolated data. */
    fclaw_debugf("Exchanging ghost patches : Setting up pointers.\n");

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

                /* Copy q and area into one contingous block */
                int pack_time_interp = time_interp && level == minlevel-1;

                double *q;
                if (gparms->manifold)
                {
                    q = (double*) e->patch_data[zz];
                    FCLAW_ASSERT(q != NULL);
                    fclaw2d_clawpatch_ghost_pack(domain,this_patch,q,
                                                 pack_time_interp);
                }
                else
                {
                    q = fclaw2d_clawpatch_get_q_timesync(domain,this_patch,
                                                         pack_time_interp);
                }
                e->patch_data[zz++] = (void*) q;
            }
        }
    }

    /* Exchange only over levels currently in use */
    if (time_interp)
    {
        int time_interp_level = minlevel-1;
        fclaw2d_domain_ghost_exchange(domain, e, time_interp_level, maxlevel);
    }
    else
    {
        fclaw2d_domain_ghost_exchange(domain, e, minlevel, maxlevel);
    }
    unpack_ghost_patches(domain,e,minlevel,maxlevel,time_interp);

    /* Count calls to this function */
    ++ddata->count_ghost_exchange;
}
