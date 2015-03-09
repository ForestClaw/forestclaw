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

#include "amr_includes.H"

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "fclaw2d_typedefs.h"
#include "fclaw2d_clawpatch.h"
#include "fclaw2d_regrid.h"

#include "fclaw2d_vtable.h"

static
void build_ghost_patches(fclaw2d_domain_t* domain)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];
        init_patch_data(ghost_patch);
        int blockno = ghost_patch->u.blockno;

        /* not clear how useful this patchno is.  In any case, it isn't
           used in defining the ClawPatch, so probably doesn't
           need to be passed in */
        int patchno = i;

        fclaw2d_clawpatch_define(domain,ghost_patch,blockno,patchno);
    }
}

/* This is called by rebuild_domain */
void fclaw2d_partition_setup(fclaw2d_domain* domain)
{
    size_t data_size =  pack_size(domain);
    fclaw2d_domain_exchange_t *e;

    /* we just created a grid by amrinit or regrid and we now need to
       allocate data to store and retrieve local boundary patches and
       remote ghost patches */
    e = fclaw2d_domain_allocate_before_exchange (domain, data_size);

    /* Store e so we can retrieve it later */
    set_domain_exchange_data(domain,e);

    /* Build patches that can be filled later with q data */
    build_ghost_patches(domain);
}

/* ------------------------------------------------------------------
   Partial exchange  - exchange_minlevel is assumed to be a time
   interpolated level.
 -------------------------------------------------------------------- */
/* Exchange_minlevel is a time interpolated level. */
void set_boundary_patch_ptrs(fclaw2d_domain_t* domain,int exchange_minlevel,
                             int exchange_maxlevel)
{
    // fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    fclaw2d_domain_exchange_t *e = get_domain_exchange_data(domain);

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

                ClawPatch *cp = get_clawpatch(this_patch);
                double *q;
                if (exchange_minlevel < level && level <= exchange_maxlevel)
                {
                    q = cp->q();
                }
                else if (level == exchange_minlevel)
                {
                    q = cp->q_time_interp();
                }
                else
                {
                    q = NULL;
                }
                e->patch_data[zz++] = (void*) q;        /* Put this patch's data location */
            }
        }
    }
}

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
            unpack_clawpatch(domain, ghost_patch,blockno, patchno, q, time_interp);
        }
    }
}

/* This is called anytime we need to update ghost patch data for certain levels
   The assumption is that the finest level is a time_interpolated level.  The
   routine 'unpack_ghost_patches' knows this, and so unpacks ghost patches to the
   correct places.
 */
void exchange_ghost_patch_data_levels(fclaw2d_domain_t* domain,
                                      int exchange_minlevel, int exchange_maxlevel)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    fclaw2d_domain_exchange_t *e = get_domain_exchange_data(domain);

    /* Do exchange to update ghost patch data */
    fclaw2d_domain_ghost_exchange(domain, e,
                                  exchange_minlevel, exchange_maxlevel);

    /* Store newly updated e->ghost_patch_data into ghost patches constructed
       locally */
    unpack_ghost_patches(domain,e, exchange_minlevel, exchange_maxlevel);

    /* Count calls to this function */
    ++ddata->count_ghost_exchange;
}

/* ------------------------------------------------------------------
   Complete exchange  - no time interpolation assumed.
 -------------------------------------------------------------------- */
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
        unpack_clawpatch(domain, ghost_patch,blockno, patchno, q, time_interp);
    }
}


/* This is called anytime we need to update ghost patch data for certain levels */
void fclaw2d_partition_exchange_all(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    fclaw2d_domain_exchange_t *e = get_domain_exchange_data(domain);

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
                ClawPatch *cp = get_clawpatch(this_patch);
                double *q = cp->q();
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



void fclaw2d_partition_domain(fclaw2d_domain_t** domain, int mode)
{
    char basename[BUFSIZ];

    // will need to access the subcyle switch
    const amr_options_t *gparms = get_domain_parms(*domain);

    // allocate memory for parallel transfor of patches
    // use data size (in bytes per patch) below.
    size_t data_size = pack_size(*domain);
    void ** patch_data = NULL;

    fclaw2d_domain_allocate_before_partition (*domain, data_size, &patch_data);

    // For all (patch i) { pack its numerical data into patch_data[i] }
    fclaw2d_domain_iterate_patches(*domain, fclaw2d_clawpatch_pack_cb,(void *) patch_data);


    // this call creates a new domain that is valid after partitioning
    // and transfers the data packed above to the new owner processors
    int exponent = gparms->subcycle && !gparms->noweightedp ? 1 : 0;
    fclaw2d_domain_t *domain_partitioned =
        fclaw2d_domain_partition (*domain, exponent);
    fclaw_bool have_new_partition = domain_partitioned != NULL;

    if (have_new_partition)
    {
        fclaw2d_domain_data_t *ddata = get_domain_data (*domain);
        fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_BUILDPATCHES]);
        fclaw2d_setup_new_domain(*domain, domain_partitioned);

	/* Stop the timer in new since, since its state is now 1. We don't care about the
	   timer in the old state.  */
        ddata = get_domain_data (domain_partitioned);
	fclaw2d_timer_stop(&ddata->timers[FCLAW2D_TIMER_BUILDPATCHES]);


        // update patch array to point to the numerical data that was received
        fclaw2d_domain_retrieve_after_partition (domain_partitioned,&patch_data);

        // TODO: for all (patch i) { unpack numerical data from patch_data[i] }
        fclaw2d_domain_iterate_patches(domain_partitioned, fclaw2d_clawpatch_unpack_cb,
                                       (void *) patch_data);

        /* then the old domain is no longer necessary */
        amrreset(domain);
        *domain = domain_partitioned;
        domain_partitioned = NULL;

        // VTK output during amrinit
        if (mode >= 0 && gparms->vtkout & 1) {
            // into timer
            fclaw2d_domain_data_t *ddata = get_domain_data (*domain);
            fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_INIT]);
            fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);

            // output
            snprintf (basename, BUFSIZ, "%s_init_level_%02d_partition",
                      gparms->prefix, mode);
            fclaw2d_output_write_vtk (*domain, basename);

            // out of timer
            fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);
            fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_INIT]);
        }

        /* internal clean up */
        fclaw2d_domain_complete(*domain);
    }

    // free the data that was used in the parallel transfer of patches
    fclaw2d_domain_free_after_partition (*domain, &patch_data);
}
