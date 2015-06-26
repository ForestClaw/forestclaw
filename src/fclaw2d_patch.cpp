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

#include "fclaw2d_patch.hpp"
#include <forestclaw2d.h>
#include <p4est_base.h>

#include <ClawPatch.H>


struct fclaw2d_patch_data
{
    ClawPatch *cp;
};


fclaw2d_patch_data_t *get_patch_data(fclaw2d_patch_t *patch)
{
    return (fclaw2d_patch_data_t *) patch->user;
}

ClawPatch* fclaw2d_patch_get_cp(fclaw2d_patch_t* this_patch)

{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    return pdata->cp;
}


void fclaw2d_patch_delete_cp(fclaw2d_patch_t* this_patch)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    if (pdata != NULL)
    {
        if (pdata->cp != NULL)
        {
            delete pdata->cp;
            pdata->cp = NULL;
        }
    }
}

/* -------------------------------------------------------------
   Put this here, since, to avoid linkage problems in fclaw2d_block.c
   ------------------------------------------------------------- */

void init_patch_data(fclaw2d_patch_t *patch)
{
    printf("init_patch_data : We shouldn't be here\n");
    exit(0);
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;
}

ClawPatch* fclaw2d_patch_new_cp(fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = new ClawPatch();
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    printf("We shouldn't be here ...\n");
    exit(0);
    pdata->cp = cp;
    return cp;
}


void fclaw2d_patch_user_data_new(fclaw2d_domain_t* domain,
                                  fclaw2d_patch_t* this_patch)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);

    /* Initialize user data */
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    this_patch->user = (void *) pdata;

    /* create new ClawPatch */
    ClawPatch *cp = new ClawPatch();
    pdata->cp = cp;
    ++ddata->count_set_clawpatch;
}

void fclaw2d_patch_user_data_delete(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t *this_patch)
{
    fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t*) this_patch->user;

    if (pdata != NULL)
    {
        fclaw2d_domain_data_t *ddata = get_domain_data(domain);
        delete pdata->cp;
        ++ddata->count_delete_clawpatch;

        FCLAW2D_FREE(pdata);
        this_patch->user = NULL;
    }
}

static
void user_data_cleanup_cb(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    fclaw2d_patch_user_data_delete(domain,this_patch);
}

void fclaw2d_patch_user_data_cleanup(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_iterate_patches(domain, user_data_cleanup_cb,
                                   (void *) NULL);
}

void init_block_and_patch_data(fclaw2d_domain_t *domain)
{
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    for (int i = 0; i < domain->num_blocks; i++)
    {
        block = &domain->blocks[i];
        init_block_data(block);
        for (int j = 0; j < block->num_patches; j++)
        {
            patch = &block->patches[j];
            patch->user = NULL;

#if 0
            fclaw2d_patch_init_user_data(domain,patch);

            init_patch_data(patch);
#endif
        }
    }
}
