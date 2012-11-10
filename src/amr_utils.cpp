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

#include "amr_utils.H"
#include "forestclaw2d.h"
#include "amr_forestclaw.H"

int pow_int(int a, int n)
{
    int b = a;
    for(int i = 0; i < n-1; i++)
    {
        b *= a;
    }
    return b;
}

// -----------------------------------------------------------------
// Initialize data
// -----------------------------------------------------------------
void init_domain_data(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = FCLAW2D_ALLOC_ZERO (fclaw2d_domain_data_t, 1);
    domain->user = (void *) ddata;
}

void init_block_data(fclaw2d_block_t *block)
{
    fclaw2d_block_data_t *bdata = FCLAW2D_ALLOC_ZERO (fclaw2d_block_data_t, 1);
    block->user = (void *) bdata;
}


void init_patch_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;
}

// -----------------------------------------------------------------
// Return pointer to user data
// -----------------------------------------------------------------
fclaw2d_domain_data_t *get_domain_data(fclaw2d_domain_t *domain)
{
    return (fclaw2d_domain_data_t *) domain->user;
}


fclaw2d_block_data_t *get_block_data(fclaw2d_block_t *block)
{
    return (fclaw2d_block_data_t *) block->user;
}


fclaw2d_patch_data_t *get_patch_data(fclaw2d_patch_t *patch)
{
    return (fclaw2d_patch_data_t *) patch->user;
}


// -----------------------------------------------------------------
// Set user data with user defined variables, etc.
// -----------------------------------------------------------------
void set_domain_data(fclaw2d_domain_t *domain, const amr_options_t *gparms)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    ddata->amropts = gparms;
}

void set_block_data(fclaw2d_block_t *block, const int mthbc[])
{
    fclaw2d_block_data_t *bdata = get_block_data(block);
    for(int i = 0; i < 4; i++)
    {
        bdata->mthbc[i] = mthbc[i];
    }
}

void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(patch);
    pdata->cp = cp;
}


// -----------------------------------------------------------------
// Some lazy helper functions that really do make things easier..
// -----------------------------------------------------------------
void allocate_user_data(fclaw2d_domain_t *domain)
{
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    init_domain_data(domain);

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


const amr_options_t* get_domain_parms(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    return ddata->amropts;
}

void set_domain_time(fclaw2d_domain_t *domain, Real time)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    ddata->curr_time = time;
}

Real get_domain_time(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    return ddata->curr_time;
}

// Will change the name of this to 'get_clawpatch' eventually
ClawPatch* get_clawpatch(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t *) patch->user;

    return pdata->cp;
}
/* end of helper functions */

const int get_refratio(fclaw2d_domain_t *domain)
{
    const amr_options_t* gparms = get_domain_parms(domain);
    return gparms->refratio;
}

// int corners_per_patch = FCLAW_CORNERS_PER_PATCH;

const int get_corners_per_patch(fclaw2d_domain_t *domain)
{
    // Number of patch corners, not the number of corners in the domain!
    return fclaw2d_domain_num_corners(domain);
}

const int get_faces_per_patch(fclaw2d_domain_t *domain)
{
    // Number of faces per patch, not the total number of faces in the domain!
    return fclaw2d_domain_num_faces(domain);
}

const int get_siblings_per_patch(fclaw2d_domain_t *domain)
{
    // Number of patch corners, not the number of corners in the domain!
    return fclaw2d_domain_num_corners(domain);
}

const int get_p4est_refineFactor(fclaw2d_domain_t *domain)
{
    return fclaw2d_domain_num_face_corners(domain);
}




static void cb_num_patches(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
  (*(int *) user)++;
}

int num_patches(fclaw2d_domain_t *domain, int level)
{
  int count = 0;
  fclaw2d_domain_iterate_level(domain, level,
                               cb_num_patches,
                               &count);
  return count;
}
