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
#include <fclaw2d_domain.h>
#include <fclaw2d_regrid.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_partition.h>

/* Also needed in amrreset */
fclaw2d_domain_exchange_t*
    fclaw2d_partition_get_exchange_data(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    return ddata->domain_exchange;
}

static
void set_exchange_data(fclaw2d_domain_t* domain,
                       fclaw2d_domain_exchange_t *e)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    ddata->domain_exchange = e;
}

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

static
void delete_ghost_patches(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];

        fclaw2d_patch_delete_cp(ghost_patch);
        fclaw2d_patch_delete_data(ghost_patch);
        ++ddata->count_delete_clawpatch;
    }
}

/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* This is called by rebuild_domain */
void fclaw2d_partition_setup(fclaw2d_domain* domain)
{
    size_t data_size =  fclaw2d_clawpatch_pack_size(domain);
    fclaw2d_domain_exchange_t *e;

    /* we just created a grid by amrinit or regrid and we now need to
       allocate data to store and retrieve local boundary patches and
       remote ghost patches */
    e = fclaw2d_domain_allocate_before_exchange (domain, data_size);

    /* Store e so we can retrieve it later */
    set_exchange_data(domain,e);

    /* Build patches that can be filled later with q data */
    build_ghost_patches(domain);
}


/* Question : Do all patches on this processor get packed? */
void fclaw2d_partition_domain(fclaw2d_domain_t** domain, int mode)
{
    char basename[BUFSIZ];

    /* will need to access the subcyle switch */
    const amr_options_t *gparms = get_domain_parms(*domain);

    /* allocate memory for parallel transfor of patches
       use data size (in bytes per patch) below. */
    size_t data_size = fclaw2d_clawpatch_pack_size(*domain);
    void ** patch_data = NULL;

    fclaw2d_domain_allocate_before_partition (*domain, data_size,
                                              &patch_data);

    /* For all (patch i) { pack its numerical data into patch_data[i] } */
    fclaw2d_domain_iterate_patches(*domain,
                                   fclaw2d_clawpatch_pack_cb,
                                   (void *) patch_data);


    /* this call creates a new domain that is valid after partitioning
       and transfers the data packed above to the new owner processors */
    int exponent = gparms->subcycle && !gparms->noweightedp ? 1 : 0;
    fclaw2d_domain_t *domain_partitioned =
        fclaw2d_domain_partition (*domain, exponent);

    fclaw_bool have_new_partition = domain_partitioned != NULL;

    if (have_new_partition)
    {
        fclaw2d_domain_data_t *ddata = get_domain_data (*domain);
        fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_BUILDPATCHES]);
        fclaw2d_regrid_new_domain_setup(*domain, domain_partitioned);

	/* Stop the timer in new since, since its state is now 1. We don't care about the
	   timer in the old state.  */
        ddata = get_domain_data (domain_partitioned);
	fclaw2d_timer_stop(&ddata->timers[FCLAW2D_TIMER_BUILDPATCHES]);


        /* update patch array to point to the numerical data that was received */
        fclaw2d_domain_retrieve_after_partition (domain_partitioned,&patch_data);

        /* TODO: for all (patch i) { unpack numerical data from patch_data[i] } */
        fclaw2d_domain_iterate_patches(domain_partitioned,
                                       fclaw2d_clawpatch_unpack_cb,
                                       (void *) patch_data);

        /* then the old domain is no longer necessary */
        amrreset(domain);
        *domain = domain_partitioned;
        domain_partitioned = NULL;

        /* VTK output during amrinit */
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

    /* free the data that was used in the parallel transfer of patches */
    fclaw2d_domain_free_after_partition (*domain, &patch_data);
}

void fclaw2d_partition_delete(fclaw2d_domain_t** domain)
{
    /* Free old parallel ghost patch data structure, must exist by construction. */
    fclaw2d_domain_data_t *ddata = get_domain_data (*domain);
    fclaw2d_domain_exchange_t *e_old = ddata->domain_exchange;

    delete_ghost_patches(*domain);
    fclaw2d_domain_free_after_exchange (*domain, e_old);
}
