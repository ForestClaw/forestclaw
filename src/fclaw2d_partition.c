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

#include <fclaw2d_partition.h>

#include <fclaw2d_convenience.h>  /* p4est domain, patch handling routines */

#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_options.h>

static
void cb_partition_pack(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *patch,
                       int blockno,
                       int patchno,
                       void *user)

{
    /* Pack everything in old domain */
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *) user;

    fclaw2d_block_t *this_block = &domain->blocks[blockno];
    int patch_num = this_block->num_patches_before + patchno;
    void* pack_data_here = (void*) ((void**)g->user)[patch_num];

    fclaw2d_patch_partition_pack(g->glob,patch,
                                 blockno,patchno,
                                 pack_data_here);
}


static
void  cb_partition_transfer(fclaw2d_domain_t * old_domain,
                            fclaw2d_patch_t * old_patch,
                            fclaw2d_domain_t * new_domain,
                            fclaw2d_patch_t * new_patch,
                            int blockno,
                            int old_patchno, int new_patchno,
                            void *user)
{
    /* Transfer data to new domain */
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *) user;
    fclaw2d_domain_data_t *ddata_old = fclaw2d_domain_get_data (old_domain);
    fclaw2d_domain_data_t *ddata_new = fclaw2d_domain_get_data (new_domain);

    if (old_patch != NULL)
    {
        FCLAW_ASSERT(old_patch->xlower == new_patch->xlower);
        FCLAW_ASSERT(old_patch->ylower == new_patch->ylower);
        FCLAW_ASSERT(old_patch->xupper == new_patch->xupper);
        FCLAW_ASSERT(old_patch->yupper == new_patch->yupper);

        new_patch->user = old_patch->user;
        old_patch->user = NULL;
        ++ddata_old->count_delete_patch;
        ++ddata_new->count_set_patch;
    }
    else
    {
        /* We need to rebuild the patch from scratch. 'user' contains
           the packed data received from remote processor. */   
        fclaw2d_domain_t *domain = new_domain;  /* get patch id in new domain */

        fclaw2d_block_t *this_block = &domain->blocks[blockno];
        int patch_num = this_block->num_patches_before + new_patchno;
        void* unpack_data_from_here = (void*) ((void**)g->user)[patch_num];

        /* pass in new_domain, since glob only contains old domain at this point
        and the both domains are needed to increment/decrement patches */
        fclaw2d_patch_partition_unpack(g->glob,new_domain,new_patch,
                                       blockno,new_patchno,unpack_data_from_here);

        /* Reason for the following two lines: the glob contains the old domain 
        which is incremented in ddata_old  but we really want to increment the 
        new domain. */
        --ddata_old->count_set_patch;
        ++ddata_new->count_set_patch;


    }
}


/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* Question : Do all patches on this processor get packed? */
void fclaw2d_partition_domain(fclaw2d_global_t* glob,
                              fclaw2d_timer_names_t running)
{
    fclaw2d_domain_t** domain = &glob->domain;
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_PARTITION]);

    /* will need to access the subcyle switch */
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    /* allocate memory for parallel transfor of patches
       use data size (in bytes per patch) below. */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_PARTITION_BUILD]);
    size_t psize = fclaw2d_patch_partition_packsize(glob);
    size_t data_size = psize;  /* Includes sizeof(data_type) */
    void ** patch_data = NULL;

    fclaw2d_domain_allocate_before_partition (*domain, data_size,
                                              &patch_data);

    /* For all (patch i) { pack its numerical data into patch_data[i] }
       Does all the data in every patch need to be copied?  */
    fclaw2d_global_iterate_patches(glob,
                                   cb_partition_pack,
                                   (void *) patch_data);
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_PARTITION_BUILD]);


    /* this call creates a new domain that is valid after partitioning
       and transfers the data packed above to the new owner processors */
    int exponent = fclaw_opt->subcycle && fclaw_opt->weighted_partition ? 1 : 0;
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_PARTITION]);
    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_stop (&glob->timers[running]);
    }
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_PARTITION_COMM]);
    fclaw2d_domain_t *domain_partitioned =
        fclaw2d_domain_partition (*domain, exponent);
    int have_new_partition = domain_partitioned != NULL;

    if (have_new_partition)
    {
        /* Do this part so we can get a pointer to the new data */
        fclaw2d_domain_setup(glob, domain_partitioned);
    }

    /* Stop the communication timer */
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_PARTITION_COMM]);
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_PARTITION]);
    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_start (&glob->timers[running]);
    }

    if (have_new_partition)
    {
        fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_PARTITION_BUILD]);

        /* update patch array to point to the numerical data that was received */
        fclaw2d_domain_retrieve_after_partition (domain_partitioned,&patch_data);

        /* New version? */
        fclaw2d_global_iterate_partitioned(glob,domain_partitioned,
                                           cb_partition_transfer,
                                           (void*) patch_data);

        /* then the old domain is no longer necessary */
        fclaw2d_domain_reset(glob);
        *domain = domain_partitioned;
        domain_partitioned = NULL;

        /* internal clean up */
        fclaw2d_domain_complete(*domain);
        fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_PARTITION_BUILD]);
    }

    /* free the data that was used in the parallel transfer of patches */
    fclaw2d_domain_free_after_partition (*domain, &patch_data);

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_PARTITION]);
}
