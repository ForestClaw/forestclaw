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

#include <forestclaw.h>
#ifndef P4_TO_P8
#include <fclaw2d_domain.h>
#include <fclaw2d_convenience.h>  /* Contains domain_destroy and others */

#include <fclaw2d_patch.h>
#include <fclaw2d_exchange.h>
#include <fclaw2d_global.h>
#else
#include <fclaw3d_domain.h>
#include <fclaw3d_convenience.h>  /* Contains domain_destroy and others */

/* when ready include <fclaw3d_patch.h> */
typedef struct fclaw3d_patch_data
{
    const fclaw3d_patch_t *real_patch;
}
fclaw3d_patch_data_t;
#endif

/* dimension-independent helper functions first */

#if 0

static fclaw2d_domain_t *
fclaw_domain_get_domain (fclaw_domain_t *d)
{
#ifndef P4_TO_P8
    return d->d.d2.domain2;
#else
    return d->d.d3.domain3;
#endif
}

#endif

void
fclaw2d_domain_iterate_cb
  (fclaw2d_domain_t * d2, fclaw2d_patch_t * patch,
   int blockno, int patchno, void *user)
{
    fclaw_domain_iterate_t *di = (fclaw_domain_iterate_t *) user;
    di->iter (di->d, (fclaw_patch_t *) patch->user, blockno, patchno,
              di->user);
}

typedef struct fcd_allocated_patch
{
#ifndef P4_TO_P8
    fclaw2d_patch_data_t pd;
#else
    fclaw3d_patch_data_t pd;
#endif
    fclaw_patch_t p;
}
fcd_allocated_patch_t;

fclaw_domain_t *
fclaw_domain_new2d (fclaw2d_domain_t * domain,
                    fclaw_domain_callback_t init, void *user)
{
    int i, j;
    fclaw2d_block_t      *block;
    fclaw2d_patch_t      *patch;
    fclaw_domain_t       *d;
    fclaw_patch_t        *p;
    fcd_allocated_patch_t *ap;

    FCLAW_ASSERT (domain != NULL && domain->mpisize > 0);

    /* allocate and set domain itself */
    d = FCLAW_ALLOC_ZERO (fclaw_domain_t, 1);
    d->dim = P4EST_DIM;
#ifndef P4_TO_P8
    d->d.d2.dmagic2 = FCLAW2D_DOMAIN_MAGIC;
    d->d.d2.domain2 = domain;
#else
    d->d.d3.dmagic3 = FCLAW3D_DOMAIN_MAGIC;
    d->d.d3.domain3 = domain;
#endif
    sc_mstamp_init (&d->pstamp, 4096 - 3 * sizeof (size_t),
                    sizeof (fcd_allocated_patch_t));

    /* iterate over all patches to initialize */
    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        for (j = 0; j < block->num_patches; ++j)
        {
            /* hook the new dimension-independent patch into storage */
            ap = (fcd_allocated_patch_t *) sc_mstamp_alloc (&d->pstamp);
            patch = block->patches + j;
            patch->user = p = &ap->p;
#ifndef P4_TO_P8
            p->pd.pd2 = &ap->pd;
            p->pd.pd2->real_patch = patch;
#else
            p->pd.pd3 = &ap->pd;
            p->pd.pd3->real_patch = patch;
#endif
            if (init != NULL) {
                init (d, p, i, j, user);
            }
        }
    }

    /* domain fully constructed */
    FCLAW_ASSERT (fclaw_domain_is_valid (d));
    return d;
}

void
fclaw_domain_destroy2d (fclaw_domain_t * d,
                        fclaw_domain_callback_t dele, void *user)
{
    int i, j;
    fclaw2d_domain_t     *domain;
    fclaw2d_block_t      *block;
    fclaw2d_patch_t      *patch;
    fclaw_patch_t        *p;

    FCLAW_ASSERT (d->dim == P4EST_DIM);
#ifndef P4_TO_P8
    domain = d->d.d2.domain2;
#else
    domain = d->d.d3.domain3;
#endif

    FCLAW_ASSERT (domain != NULL && domain->mpisize > 0);
    FCLAW_ASSERT (d->du.count_set_patch == d->du.count_delete_patch);

    /* iterate over all patches to deinitialize */
    if (dele != NULL) {
        for (i = 0; i < domain->num_blocks; ++i)
        {
            block = domain->blocks + i;
            for (j = 0; j < block->num_patches; ++j)
            {
                /* free the new dimension-independent patch into storage */
                patch = block->patches + j;
                p = (fclaw_patch_t *) patch->user;
                dele (d, p, i, j, user);
            }
        }
    }
    fclaw2d_domain_destroy (domain);
    sc_mstamp_reset (&d->pstamp);
    FCLAW_FREE (d);
}

#ifndef P4_TO_P8

/* we're holding back with 3d counterparts
   since much of this will move into fclaw_domain.c */
/* below follows the previous code unchanged */

void fclaw2d_domain_data_new(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t* ddata = (fclaw2d_domain_data_t*) domain->user;
    ddata = FCLAW_ALLOC_ZERO (fclaw2d_domain_data_t, 1);
    domain->user = ddata;

    ddata->count_set_patch = ddata->count_delete_patch = 0;
    
    ddata->domain_exchange = NULL;
    ddata->domain_indirect = NULL;
}

void fclaw2d_domain_data_delete(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t* ddata = (fclaw2d_domain_data_t*) domain->user;

    FCLAW_FREE (ddata);
    domain->user = NULL;
}

fclaw2d_domain_data_t *fclaw2d_domain_get_data(fclaw2d_domain_t *domain)
{
    return (fclaw2d_domain_data_t *) domain->user;
}


void fclaw2d_domain_setup(fclaw2d_global_t* glob,
                          fclaw2d_domain_t* new_domain)
{
    fclaw2d_domain_t *old_domain = glob->domain;
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
        fclaw2d_domain_data_new(new_domain);
    }
    fclaw_global_infof("Done\n");
}


void fclaw2d_domain_reset(fclaw2d_global_t* glob)
{
    fclaw2d_domain_t** domain = &glob->domain;
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (*domain);
    int i, j;

    for(i = 0; i < (*domain)->num_blocks; i++)
    {
        fclaw2d_block_t *block = (*domain)->blocks + i;

        for(j = 0; j < block->num_patches; j++)
        {
            /* This is here to delete any patches created during
               initialization, and not through regridding */
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_data_delete(glob,patch);
        }
        block->user = NULL;

    }

    if (ddata->domain_exchange != NULL)
    {
        fclaw2d_exchange_delete(glob);
    }

    /* Output memory discrepancy for the ClawPatch */
    if (ddata->count_set_patch != ddata->count_delete_patch)
    {
        printf ("[%d] This domain had Clawpatch set %d and deleted %d times\n",
                (*domain)->mpirank,
                ddata->count_set_patch, ddata->count_delete_patch);
    }

    fclaw2d_domain_data_delete(*domain);  // Delete allocated pointers to set of functions.

    fclaw2d_domain_destroy(*domain);
    *domain = NULL;
}


void fclaw2d_domain_iterate_level_mthread (fclaw2d_domain_t * domain, int level,
                                           fclaw2d_patch_callback_t pcb, void *user)
{
#if (_OPENMP)
    int i, j;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

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
    fclaw_global_essentialf("fclaw2d_patch_iterator_mthread : We should not be here\n");
#endif
}

#endif /* !P4_TO_P8 */
