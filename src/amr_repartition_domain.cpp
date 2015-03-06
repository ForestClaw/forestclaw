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
#include "fclaw2d_regrid.H"

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

        set_clawpatch(domain,ghost_patch,blockno,patchno);
    }
}

void delete_ghost_patches(fclaw2d_domain_t* domain)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];

        delete_clawpatch(domain, ghost_patch);
        delete_patch_data(ghost_patch);
    }
}

/* This is called by rebuild_domain */
static
void setup_parallel_ghost_exchange(fclaw2d_domain_t* domain)
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

#if 0
/* This is called anytime we need to update ghost patch data */
void exchange_ghost_patch_data(fclaw2d_domain_t* domain, fclaw_bool time_interp)
{
    exchange_ghost_patch_data_levels (domain, time_interp,
                                      domain->global_minlevel, domain->global_maxlevel);
}
#endif


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
void exchange_ghost_patch_data_all(fclaw2d_domain_t* domain)
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



/* ------------------------------------------------------------------
   Repartition and rebuild new domains, or construct initial domain
 -------------------------------------------------------------------- */

/* Build initial set of patches */
static
void cb_build_patches(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    set_clawpatch(domain,this_patch,this_block_idx,this_patch_idx);

    if (vt.patch_setup != NULL)
    {
        vt.patch_setup(domain,this_patch,this_block_idx,this_patch_idx);
    }
}


static
void cb_pack_patches(fclaw2d_domain_t *domain,
                     fclaw2d_patch_t *this_patch,
                     int this_block_idx,
                     int this_patch_idx,
                     void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    pack_clawpatch(this_patch,patch_data);
}

static
void cb_unpack_patches(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    fclaw_bool time_interp = fclaw_false;
    unpack_clawpatch(domain,this_patch,this_block_idx,this_patch_idx,
                     patch_data,time_interp);
}



void rebuild_domain(fclaw2d_domain_t* old_domain, fclaw2d_domain_t* new_domain)
{
    fclaw_global_infof("Rebuilding domain\n");

    const amr_options_t *gparms = get_domain_parms(old_domain);
    double t = get_domain_time(old_domain);

    // Allocate memory for user data types (but they don't get set)
    init_domain_data(new_domain);

    copy_domain_data(old_domain,new_domain);

    // Why isn't this done in copy_domain_data?
    set_domain_time(new_domain,t);

    // Allocate memory for new blocks and patches.
    init_block_and_patch_data(new_domain);

    // Physical BCs are needed in boundary level exchange
    // Assume only one block, since we are assuming mthbc
    int num = new_domain->num_blocks;
    for (int i = 0; i < num; i++)
    {
        fclaw2d_block_t *block = &new_domain->blocks[i];
        // This is kind of dumb for now, since block won't in general
        // have the same physical boundary conditions types.
        set_block_data(block,gparms->mthbc);
    }

    fclaw_global_infof("  -- Rebuilding patches ... ");

    fclaw2d_domain_iterate_patches(new_domain, cb_build_patches,(void *) NULL);

    fclaw_global_infof("Done\n");

    // Set up the parallel ghost patch data structure.
    fclaw_global_infof("  -- Setting up parallel ghost exchange ... ");

    setup_parallel_ghost_exchange(new_domain);

    fclaw_global_infof("Done\n");
}

void build_initial_domain(fclaw2d_domain_t* domain)
{
    fclaw_global_infof("Building initial domain\n");

    const amr_options_t *gparms = get_domain_parms(domain);

    // init_domain_data(*domain) is not called here, because it is
    // called in <main>.cpp.  This allows the user to grab gparms,
    // setup_problem(), etc before getting here .

    // Allocate memory for new blocks and patches.
    init_block_and_patch_data(domain);

    // Physical BCs are needed in boundary level exchange
    // Assume only one block, since we are assuming mthbc
    int num = domain->num_blocks;
    for (int i = 0; i < num; i++)
    {
        // This assumes that each block has the same physical
        // boundary conditions, which doesn't make much sense...
        fclaw2d_block_t *block = &domain->blocks[i];
        set_block_data(block,gparms->mthbc);
    }

    // Construct new patches
    fclaw_global_infof("  -- Rebuilding patches ... ");

    fclaw2d_domain_iterate_patches(domain, cb_build_patches,(void *) NULL);

    fclaw_global_infof("Done\n");

    // Set up the parallel ghost patch data structure.
    fclaw_global_infof("  -- Setting up parallel ghost exchange ... ");

    setup_parallel_ghost_exchange(domain);

    fclaw_global_infof("Done\n");

}


void repartition_domain(fclaw2d_domain_t** domain, int mode)
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
    fclaw2d_domain_iterate_patches(*domain, cb_pack_patches,(void *) patch_data);


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
        rebuild_domain(*domain, domain_partitioned);

	/* Stop the timer in new since, since its state is now 1. We don't care about the
	   timer in the old state.  */
        ddata = get_domain_data (domain_partitioned);
	fclaw2d_timer_stop(&ddata->timers[FCLAW2D_TIMER_BUILDPATCHES]);


        // update patch array to point to the numerical data that was received
        fclaw2d_domain_retrieve_after_partition (domain_partitioned,&patch_data);

        // TODO: for all (patch i) { unpack numerical data from patch_data[i] }
        fclaw2d_domain_iterate_patches(domain_partitioned, cb_unpack_patches,
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


/* ------------------------------------------------------------------
   Print out diagnostic information
 -------------------------------------------------------------------- */

static
void cb_proc_info (fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    int level = this_patch->level;
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
        domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);

    fclaw_debugf("%5d %5d %5d\n",(int) patch_num,level,this_patch_idx);
}



void amr_print_patches_and_procs(fclaw2d_domain_t *domain)
{
        fclaw2d_domain_iterate_patches(domain, cb_proc_info,(void *) NULL);
}
