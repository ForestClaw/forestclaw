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

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "fclaw2d_typedefs.h"
#include "amr_regrid.H"

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

static
void unpack_ghost_patches(fclaw2d_domain_t* domain, fclaw2d_domain_exchange_t *e)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];
        int blockno = ghost_patch->u.blockno;

        /* not clear how useful this patchno is.  In any case, it isn't
           used in defining the ClawPatch, so probably doesn't
           need to be passed in */
        int patchno = i;

        /* access data stored on remote procs */
        double *q = (double*) e->ghost_data[patchno];

        unpack_clawpatch(domain, ghost_patch,blockno, patchno, q);

        // should not be needed here
        // set_clawpatch(domain,ghost_patch,blockno,patchno);
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



/* This is called anytime we need to update ghost patch data */
void exchange_ghost_patch_data(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_exchange_t *e = get_domain_exchange_data(domain);

    /* Store local boundary data */
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
                e->patch_data[zz++] = (void*) q;        /* Put this patch's data location */
            }
        }
    }

    /* Do exchange to update ghost patch data */
    fclaw2d_domain_ghost_exchange(domain, e);

    /* Store newly updated e->ghost_patch_data into ghost patches constructed
       locally */
    unpack_ghost_patches(domain,e);
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
    // Grid doesn't change
    set_clawpatch(domain,this_patch,this_block_idx,this_patch_idx);

    // Setup new patch using solver specific routine
    fclaw2d_solver_functions_t *sf = get_solver_functions(domain);
    (sf->f_patch_setup)(domain,this_patch,this_block_idx,this_patch_idx);
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

    unpack_clawpatch(domain,this_patch,this_block_idx,this_patch_idx,
                     patch_data);
}



void rebuild_domain(fclaw2d_domain_t* old_domain, fclaw2d_domain_t* new_domain)
{
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

    fclaw2d_domain_iterate_patches(new_domain, cb_build_patches,(void *) NULL);

    // Set up the parallel ghost patch data structure.
    fclaw2d_domain_exchange_t *e_old = get_domain_exchange_data(new_domain);
    if (e_old != NULL)
    {
        fclaw2d_domain_free_after_exchange (new_domain, e_old);
    }

    setup_parallel_ghost_exchange(new_domain);
}

void build_initial_domain(fclaw2d_domain_t* domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);

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
    fclaw2d_domain_iterate_patches(domain, cb_build_patches,(void *) NULL);

    // Set up the parallel ghost patch data structure.
    fclaw2d_domain_exchange_t *e_old = get_domain_exchange_data(domain);
    if (e_old != NULL)
    {
        fclaw2d_domain_free_after_exchange (domain, e_old);
    }
    setup_parallel_ghost_exchange(domain);
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
    fclaw2d_domain_t *domain_partitioned =
        fclaw2d_domain_partition (*domain, gparms->subcycle ? 1 : 0);
    fclaw_bool have_new_partition = domain_partitioned != NULL;

    if (have_new_partition)
    {
        rebuild_domain(*domain, domain_partitioned);

        // update patch array to point to the numerical data that was received
        fclaw2d_domain_retrieve_after_partition (domain_partitioned,&patch_data);

        // TODO: for all (patch i) { unpack numerical data from patch_data[i] }
        fclaw2d_domain_iterate_patches(domain_partitioned, cb_unpack_patches,
                                       (void *) patch_data);

        /* then the old domain is no longer necessary */
        amrreset(domain);
        *domain = domain_partitioned;

        // VTK output during amrinit
        if (mode >= 0 && gparms->vtkout & 1) {
            snprintf (basename, BUFSIZ, "init_level_%02d_partition", mode);
            amr_output_write_vtk (*domain, basename);
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
    int mpirank = domain->mpirank;
    int level = this_patch->level;
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
        domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);
    printf("%5d %5d %5d %5d\n",mpirank, (int) patch_num,level,this_patch_idx);
}



void amr_print_patches_and_procs(fclaw2d_domain_t *domain)
{
        fclaw2d_domain_iterate_patches(domain, cb_proc_info,(void *) NULL);
}
