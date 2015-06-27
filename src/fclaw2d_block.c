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

#include <fclaw2d_global.h>
#include <fclaw2d_forestclaw.h>

fclaw2d_block_data_t *fclaw2d_block_get_data(fclaw2d_block_t *block)
{
    return (fclaw2d_block_data_t *) block->user;
}

void fclaw2d_block_set_data(fclaw2d_block_t *block, const int mthbc[])
{
    fclaw2d_block_data_t *bdata = fclaw2d_block_get_data(block);
    int i;

    for (i = 0; i < 4; i++)
    {
        bdata->mthbc[i] = mthbc[i];
    }
}

void fclaw2d_block_get_block_boundary(fclaw2d_domain_t * domain,
                                      fclaw2d_patch_t * patch,
                                      fclaw_bool *intersects_block)
{
    int iside;

    for (iside = 0; iside < NumFaces; iside++)
    {
        int iface_flags = fclaw2d_patch_block_face_flags[iside];
        int is_block_face = (patch->flags & iface_flags) != 0;

        /* True for physical and block boundaries across a face */
        intersects_block[iside] = is_block_face;
    }
}

static
void block_data_new(fclaw2d_block_t *block)
{
    fclaw2d_block_data_t *bdata = FCLAW2D_ALLOC_ZERO (fclaw2d_block_data_t, 1);
    block->user = (void *) bdata;
}

void fclaw2d_block_data_new(fclaw2d_domain_t *domain)
{
    fclaw2d_block_t *block;
    int i;

    for (i = 0; i < domain->num_blocks; i++)
    {
        block = &domain->blocks[i];
        block_data_new(block);
    }
}
