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

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static
void  cb_partition_transfer(fclaw2d_domain_t * old_domain,
                            fclaw2d_patch_t * old_patch,
                            fclaw2d_domain_t * new_domain,
                            fclaw2d_patch_t * new_patch,
                            int blockno,
                            int old_patchno, int new_patchno,
                            void *user)
{
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
        ++ddata_old->count_delete_clawpatch;
        ++ddata_new->count_set_clawpatch;
    }
    else
    {
        /* We need to rebuild the patch from scratch. 'user' contains
           the packed data received from remote processor. */
        cb_fclaw2d_clawpatch_partition_unpack(new_domain,new_patch,
                                              blockno,new_patchno,user);
    }
}


/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* Question : Do all patches on this processor get packed? */
void fclaw2d_partition_domain(fclaw2d_domain_t** domain, int mode)
{
    char basename[BUFSIZ];

    /* will need to access the subcyle switch */
    const amr_options_t *gparms = get_domain_parms(*domain);

    /* allocate memory for parallel transfor of patches
       use data size (in bytes per patch) below. */
    size_t data_size = fclaw2d_clawpatch_partition_packsize(*domain);
    void ** patch_data = NULL;

    fclaw2d_domain_allocate_before_partition (*domain, data_size,
                                              &patch_data);

    /* For all (patch i) { pack its numerical data into patch_data[i] }
       Does all the data in every patch need to be copied?  */
    fclaw2d_domain_iterate_patches(*domain,
                                   cb_fclaw2d_clawpatch_partition_pack,
                                   (void *) patch_data);


    /* this call creates a new domain that is valid after partitioning
       and transfers the data packed above to the new owner processors */
    int exponent = gparms->subcycle && !gparms->noweightedp ? 1 : 0;
    fclaw2d_domain_t *domain_partitioned =
        fclaw2d_domain_partition (*domain, exponent);

    fclaw_bool have_new_partition = domain_partitioned != NULL;

    if (have_new_partition)
    {
        fclaw2d_domain_setup(*domain, domain_partitioned);

        /* update patch array to point to the numerical data that was received */
        fclaw2d_domain_retrieve_after_partition (domain_partitioned,&patch_data);

        /* TODO: for all (patch i) { unpack numerical data from patch_data[i] } */
        fclaw2d_domain_data_t* ddata = fclaw2d_domain_get_data(domain_partitioned);
        fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_BUILDPARTITION]);

#if 0
        /* Old version */
        fclaw2d_domain_iterate_patches(domain_partitioned,
                                       cb_fclaw2d_clawpatch_partition_unpack,
                                       (void *) patch_data);
#endif

        /* New version? */
        fclaw2d_domain_iterate_partitioned(*domain,domain_partitioned,
                                           cb_partition_transfer,
                                           (void*) patch_data);

        fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_BUILDPARTITION]);

        /* then the old domain is no longer necessary */
        fclaw2d_domain_reset(domain);
        *domain = domain_partitioned;
        domain_partitioned = NULL;

        /* VTK output during amrinit */
        if (mode >= 0 && gparms->vtkout & 1) {
            // into timer
            fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (*domain);
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

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
