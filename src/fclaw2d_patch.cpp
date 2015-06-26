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


void init_patch_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;
}

void fclaw2d_patch_delete_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pd = (fclaw2d_patch_data_t*) patch->user;
    FCLAW2D_FREE(pd);
}

fclaw2d_patch_data_t *get_patch_data(fclaw2d_patch_t *patch)
{
    return (fclaw2d_patch_data_t *) patch->user;
}

ClawPatch* fclaw2d_patch_get_cp(fclaw2d_patch_t* this_patch)

{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    return pdata->cp;
}

ClawPatch* fclaw2d_patch_new_cp(fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = new ClawPatch();
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    pdata->cp = cp;
    return cp;
}

void fclaw2d_patch_delete_cp(fclaw2d_patch_t* this_patch)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    FCLAW_ASSERT(pdata->cp != NULL);
    delete pdata->cp;
    pdata->cp = NULL;
}

/* -------------------------------------------------------------
   Put this here, since, to avoid linkage problems in fclaw2d_block.c
   ------------------------------------------------------------- */
void init_block_and_patch_data(fclaw2d_domain_t *domain)
{
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    // init_domain_data(domain);

    for (int i = 0; i < domain->num_blocks; i++)
    {
        block = &domain->blocks[i];
        init_block_data(block);
        for (int j = 0; j < block->num_patches; j++)
        {
            patch = &block->patches[j];
            init_patch_data(patch);
        }
    }
}



#if 0
static void cb_num_patches(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
  (*(int *) user)++;
}

int num_patches(fclaw2d_domain_t *domain, int level, int include_shadow)
{
    int count = 0;
    if (include_shadow == 0)
    {
        fclaw2d_domain_iterate_level(domain, level,
                                     cb_num_patches,
                                     &count);
    }
    else
    {
        // Include shadow patches
    }
    return count;
}
#endif
