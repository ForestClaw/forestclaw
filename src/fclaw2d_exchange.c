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

#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw_options.h>

#ifndef P4_TO_P8

#include <fclaw2d_exchange.h>

#include <fclaw2d_patch.h>


#else

#include <fclaw3d_exchange.h>

#include <fclaw3d_patch.h>

#include <fclaw2d_to_3d.h>

#define d2 d3

#endif

/* Also needed in fclaw2d_domain_reset */
/* edit: right now this is not the case, switching to static */
static
fclaw2d_domain_exchange_t*
    get_exchange_data(fclaw_global_t* glob)
{
    fclaw2d_domain_data_t *ddata = glob->domain->d2;
    return ddata->domain_exchange;
}

/* Should these be access functions in domain?  */
static
void set_exchange_data(fclaw_global_t* glob,
                       fclaw2d_domain_exchange_t *e)
{
    fclaw2d_domain_data_t *ddata = glob->domain->d2;
    ddata->domain_exchange = e;
}

static
fclaw2d_domain_indirect_t*
    get_indirect_data(fclaw_global_t* glob)
{
    fclaw2d_domain_data_t *ddata = glob->domain->d2;
    return ddata->domain_indirect;
}

static
void set_indirect_data(fclaw_global_t* glob,
                       fclaw2d_domain_indirect_t *ind)
{
    fclaw2d_domain_data_t *ddata = glob->domain->d2;
    ddata->domain_indirect = ind;
}

static
void build_remote_ghost_patches(fclaw_global_t* glob)
{
    fclaw_domain_t *domain = glob->domain;
    const fclaw_options_t *gparms = fclaw_get_options(glob);

    fclaw_infof("[%d] Number of ghost patches : %d\n",
                            domain->mpirank,domain->num_ghost_patches);
    int blockno, patchno;
    fclaw_patch_t *ghost_patch;
    fclaw_build_mode_t build_mode;
    if (gparms->ghost_patch_pack_area)
    {
        build_mode = FCLAW_BUILD_FOR_GHOST_AREA_PACKED;
    }
    else
    {
        build_mode = FCLAW_BUILD_FOR_GHOST_AREA_COMPUTED;
    }

    int i;
    for(i = 0; i < domain->num_ghost_patches; i++)
    {
        ghost_patch = &domain->ghost_patches[i];

        blockno = ghost_patch->u.blockno;

        /* not clear how useful this patchno is.  In any case, it isn't
           used in defining the ClawPatch, so probably doesn't
           need to be passed in */
        patchno = i;

        fclaw2d_patch_remote_ghost_build(glob,ghost_patch,blockno,
                                         patchno, build_mode);
    }
    fclaw_infof("[%d] Done building remote ghost patches : %d\n",
                            domain->mpirank,domain->num_ghost_patches);
}

static
void delete_remote_ghost_patches(fclaw_global_t* glob)
{
    fclaw_domain_t *domain = glob->domain;
    int i;
    for(i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw_patch_t* ghost_patch = &domain->ghost_patches[i];
        
        fclaw2d_patch_remote_ghost_delete(glob,ghost_patch);
    }
}


static void
unpack_remote_ghost_patches(fclaw_global_t* glob,
                            fclaw2d_domain_exchange_t *e,
                            int minlevel,
                            int maxlevel,
                            int time_interp)
{
    fclaw_domain_t *domain = glob->domain;

    int i;
    for(i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw_patch_t* ghost_patch = &domain->ghost_patches[i];
        int level = ghost_patch->level;

        if (level >= minlevel-1)
        {
            int blockno = ghost_patch->u.blockno;

            int patchno = i;

            /* access data stored on remote procs. */
            void *q = e->ghost_data[patchno];

            int unpack_to_timeinterp_patch=0;
            if (time_interp && level == minlevel-1)
            {
                unpack_to_timeinterp_patch = 1;
            }
            fclaw2d_patch_remote_ghost_unpack(glob, ghost_patch, blockno,
                                              patchno, q, unpack_to_timeinterp_patch);
        }
    }
}

/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* This is called whenever a new domain is created (initialize, regrid) */
void fclaw2d_exchange_setup(fclaw_global_t* glob,
                            fclaw_timer_names_t running)
{
    fclaw_domain_t* domain = glob->domain;
    /* Time spend in build here is negligible and is included in regrid */
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);

    size_t psize = fclaw2d_patch_ghost_packsize(glob);
    size_t data_size = psize;  /* Includes sizeof(datatype), i.e. sizeof(double) */
    fclaw2d_domain_exchange_t *e;

    /* we just created a grid by fclaw2d_initialize or fclaw2d_regrid
       and we now need to allocate data to store and retrieve local
       boundary patches and remote ghost patches */
    e = fclaw2d_domain_allocate_before_exchange (domain, data_size);

    /* Store e so we can retrieve it later */
    set_exchange_data(glob,e);

    /* Store locations of on-proc boundary patches that will be communicated
       to neighboring remote procs.  The pointer stored in e->patch_data[]
       will point to either the patch data itself, or to a newly allocated
       contiguous memory block that stores qdata and area (for computations
       on manifolds).*/
    int zz = 0;
    int nb;
    for (nb = 0; nb < domain->num_blocks; ++nb)
    {
        int np;
        for (np = 0; np < domain->blocks[nb].num_patches; ++np)
        {
            if (domain->blocks[nb].patches[np].flags &
                FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY)
            {                
                fclaw2d_patch_local_ghost_alloc(glob, &e->patch_data[zz++]);
            }
        }
    }

    /* ---------------------------------------------------------
       Start send for ghost patch meta-data needed to handle
       multi-proc corner case

       Note that this is a parallel communication, but only needs
       to happen once when setting up the exchange.
       ------------------------------------------------------- */
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);
    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_stop (&glob->timers[running]);
    }

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);
    fclaw2d_domain_indirect_t *ind =
        fclaw2d_domain_indirect_begin(domain);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);

    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_start (&glob->timers[running]);
    }

    /* Do some some work that we hope to hide by communication above.  */
    set_indirect_data(glob,ind);

    /* Build ghost patches from neighboring remote processors.  These will be
       filled later with q data and the area, if we are on a manifold */

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);
    build_remote_ghost_patches(glob);
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);

    /* ---------------------------------------------------------
       Receive ghost patch meta data from send initiated above.
       ------------------------------------------------------- */
    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_stop (&glob->timers[running]);
    }

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);
    fclaw2d_domain_indirect_end(domain,ind);
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);

    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_start (&glob->timers[running]);
    }
}

void fclaw2d_exchange_delete(fclaw_global_t* glob)
{
    fclaw_domain_t** domain = &glob->domain;
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);

    /* Free old parallel ghost patch data structure, must exist by construction. */
    fclaw2d_domain_exchange_t *e_old = get_exchange_data(glob);

    /* Free contiguous memory, if allocated.  If no memory was allocated
       (because we are not storing the area), nothing de-allocated. */
    if (e_old != NULL)
    {
        /* Delete local patches which are passed to other procs */
        int zz = 0;
        int nb;
        for (nb = 0; nb < (*domain)->num_blocks; ++nb)
        {
            int np;
            for (np = 0; np < (*domain)->blocks[nb].num_patches; ++np)
            {
                if ((*domain)->blocks[nb].patches[np].flags &
                    FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY)
                {
                    fclaw2d_patch_local_ghost_free(glob,&e_old->patch_data[zz++]);
                }
            }
        }
    }

    /* Delete ghost patches from remote neighboring patches */
    delete_remote_ghost_patches(glob);
    fclaw2d_domain_free_after_exchange (*domain, e_old);

    /* Destroy indirect data needed to communicate between ghost patches
       from different procs */
    fclaw2d_domain_indirect_t* ind_old = get_indirect_data(glob);
    fclaw2d_domain_indirect_destroy(*domain,ind_old);
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);
}


/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */

/* This is called whenever all time levels are time synchronized. */
void fclaw2d_exchange_ghost_patches_begin(fclaw_global_t* glob,
                                          int minlevel,
                                          int maxlevel,
                                          int time_interp,
                                          fclaw_timer_names_t running)
{
    fclaw_domain_t* domain = glob->domain;
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);

    fclaw2d_domain_exchange_t *e = get_exchange_data(glob);

    /* Pack local data into on-proc patches at the parallel boundary that
       will be shipped of to other processors. */
    int zz = 0;
    int nb;
    for (nb = 0; nb < domain->num_blocks; ++nb)
    {
        int np;
        for (np = 0; np < domain->blocks[nb].num_patches; ++np)
        {
            if (domain->blocks[nb].patches[np].flags &
                FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY)
            {
                fclaw_patch_t *this_patch = &domain->blocks[nb].patches[np];
                int level = this_patch->level;
                FCLAW_ASSERT(level <= maxlevel);

                int pack_time_interp = time_interp && level == minlevel-1;

                /* Pack q and area into one contingous block */
                fclaw2d_patch_local_ghost_pack(glob,this_patch,
                                               e->patch_data[zz++],
                                               pack_time_interp);
            }
        }
    }

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);

    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_stop (&glob->timers[running]);
    }
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);
    /* Exchange only over levels currently in use */
    if (time_interp)
    {
        int time_interp_level = minlevel-1;

        fclaw2d_domain_ghost_exchange_begin
            (domain, e,time_interp_level, maxlevel);
    }
    else
    {
        fclaw2d_domain_ghost_exchange_begin(domain, e, minlevel, maxlevel);

    }
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);
    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_start (&glob->timers[running]);
    }
}

/* This is called whenever all time levels are time synchronized. */
void fclaw2d_exchange_ghost_patches_end(fclaw_global_t* glob,
                                        int minlevel,
                                        int maxlevel,
                                        int time_interp,
                                        fclaw_timer_names_t running)
{
    fclaw_domain_t* domain = glob->domain;
    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_stop (&glob->timers[running]);
    }
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);

    fclaw2d_domain_exchange_t *e = get_exchange_data(glob);

    /* Exchange only over levels currently in use */
    if (time_interp)
    {
        fclaw2d_domain_ghost_exchange_end (domain, e);
    }
    else
    {
        fclaw2d_domain_ghost_exchange_end (domain, e);
    }

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_COMM]);
    if (running != FCLAW_TIMER_NONE)
    {
        fclaw_timer_start (&glob->timers[running]);
    }


    fclaw_timer_start (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);
    /* Unpack data from remote patches to corresponding ghost patches
       stored locally */
    unpack_remote_ghost_patches(glob,e,minlevel,maxlevel,time_interp);

    /* Count calls to this function */
    ++glob->count_ghost_exchange;
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_GHOSTPATCH_BUILD]);
}
