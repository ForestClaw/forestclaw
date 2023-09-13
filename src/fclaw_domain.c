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

#include <fclaw_convenience.h>  /* Contains domain_destroy and others */
#include <fclaw_patch.h>
#include <fclaw_exchange.h>

#include <fclaw_patch.h>


/* we're holding back with 3d counterparts
   since much of this will move into fclaw_domain.c */
/* below follows the previous code unchanged */
/* edit: do it anyway for now to support p8/swirl */

void fclaw_domain_setup(fclaw_global_t* glob,
                          fclaw_domain_t* new_domain)
{
    fclaw_domain_t *old_domain = glob->domain;
    double t;

    if (old_domain == new_domain)
    {
        fclaw_global_infof("Building initial domain\n");
        t = 0;
        glob->curr_time = t;//new_domain        
    }
    else
    {
        fclaw_global_infof("Rebuilding  domain\n");
    }
    fclaw_global_infof("Done\n");
}

void fclaw_domain_reset(fclaw_global_t* glob)
{
    fclaw_domain_t** domain = &glob->domain;
    int i, j;

    for(i = 0; i < (*domain)->num_blocks; i++)
    {
        fclaw_block_t *block = (*domain)->blocks + i;

        for(j = 0; j < block->num_patches; j++)
        {
            /* This is here to delete any patches created during
               initialization, and not through regridding */
            fclaw_patch_t *patch = block->patches + j;
            fclaw_patch_data_delete(glob,patch);
        }
        block->user = NULL;
    }

    // TODO provide function in fclaw_exchange
    if (((*domain)->dim == 2 && (*domain)->d2->exchange != NULL)
        || ((*domain)->dim == 3 && (*domain)->d3->exchange != NULL))
    {
        fclaw_exchange_delete(glob);
    }

    FCLAW_ASSERT((*domain)->count_set_patch == (*domain)->count_delete_patch);

    fclaw_domain_destroy(*domain);
    *domain = NULL;
}

void fclaw_domain_iterate_level_mthread (fclaw_domain_t * domain, int level,
                                           fclaw_patch_callback_t pcb, void *user)
{
#if (_OPENMP)
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

    for (i = 0; i < domain->num_blocks; i++)
    {
        block = domain->blocks + i;
#pragma omp parallel for private(patch,j)
        for (j = 0; j < block->num_patches; j++)
        {
            patch = block->patches + j;
            if (patch->level == level)
            {
                pcb (domain, patch, i, j, user);
            }
        }
    }
#else
#ifndef P4_TO_P8
    fclaw_global_essentialf("fclaw2d_patch_iterator_mthread: We should not be here\n");
#else
    fclaw_global_essentialf("fclaw3d_patch_iterator_mthread: We should not be here\n");
#endif
#endif
}
