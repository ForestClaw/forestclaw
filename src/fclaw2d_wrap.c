/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw2d_wrap.h>
#include <fclaw2d_defs.h>

#else

#include <fclaw3d_wrap.h>
#include <fclaw3d_defs.h>

#endif

static
void copy_patch(fclaw_patch_t* patch, fclaw2d_patch_t* patch_2d)
{
    patch->refine_dim = FCLAW2D_SPACEDIM;
#ifndef P4_TO_P8
    patch->d2 = patch_2d;
#else
    patch->d3 = patch_2d;
#endif
    patch->xlower = patch_2d->xlower;
    patch->xupper = patch_2d->xupper;
    patch->ylower = patch_2d->ylower;
    patch->yupper = patch_2d->yupper;
#ifndef P4_TO_P8
    patch->zlower = 0;
    patch->zupper = 1;
#else
    patch->zlower = patch_2d->zlower;
    patch->zupper = patch_2d->zupper;
#endif


    patch->level = patch_2d->level;
    patch_2d->user = patch;
}

static
void copy_block(fclaw_block_t* block, fclaw2d_block_t* block_2d)
{
    block->refine_dim = FCLAW2D_SPACEDIM;
#ifndef P4_TO_P8
    block->d2 = block_2d;
#else
    block->d3 = block_2d;
#endif

    block->num_patches = block_2d->num_patches;
    block->num_patches_before = block_2d->num_patches_before;

    block->vertices = block_2d->vertices;

    block->patches = FCLAW_ALLOC_ZERO(fclaw_patch_t, block->num_patches);
    for(int patchno = 0; patchno < block->num_patches; patchno++)
    {
        copy_patch(&block->patches[patchno], &block_2d->patches[patchno]);
    }
}

fclaw_domain_t* fclaw_domain_wrap_2d(fclaw2d_domain_t* domain_2d)
{
    if(domain_2d == NULL)
    {
        return NULL;
    }

    fclaw_domain_t* domain = FCLAW_ALLOC_ZERO(fclaw_domain_t, 1);
    domain->refine_dim = FCLAW2D_SPACEDIM;
#ifndef P4_TO_P8
    domain->d2 = domain_2d;
#else
    domain->d3 = domain_2d;
#endif

    domain->mpicomm = domain_2d->mpicomm;
    domain->mpisize = domain_2d->mpisize;
    domain->mpirank = domain_2d->mpirank;

    domain->local_num_patches = domain_2d->local_num_patches;
    domain->local_max_patches = domain_2d->local_max_patches;
    domain->local_minlevel = domain_2d->local_minlevel;
    domain->local_maxlevel = domain_2d->local_maxlevel;
    domain->global_num_patches = domain_2d->global_num_patches;
    domain->global_num_patches_before = domain_2d->global_num_patches_before;
    domain->global_minlevel = domain_2d->global_minlevel;
    domain->global_maxlevel = domain_2d->global_maxlevel;

    domain->num_blocks = domain_2d->num_blocks;

    domain->blocks = FCLAW_ALLOC_ZERO(fclaw_block_t, domain->num_blocks);

    for(int blockno = 0; blockno < domain->num_blocks; ++blockno)
    {
        copy_block(&domain->blocks[blockno], &domain_2d->blocks[blockno]);
    }

    domain->num_exchange_patches = domain_2d->num_exchange_patches;

    domain->num_ghost_patches = domain_2d->num_ghost_patches;
    domain->ghost_patches = FCLAW_ALLOC_ZERO(fclaw_patch_t, domain->num_ghost_patches);
    for(int patchno = 0; patchno < domain->num_ghost_patches; patchno++)
    {
        copy_patch(&domain->ghost_patches[patchno], &domain_2d->ghost_patches[patchno]);
    }

    domain->attributes = domain_2d->attributes;

    domain_2d->user = domain;

    return domain;
}

static
fclaw_domain_t* get_domain(fclaw2d_domain_t* domain_2d)
{
    FCLAW_ASSERT(domain_2d->user != NULL);
    return (fclaw_domain_t*) domain_2d->user;
}

static
fclaw_patch_t* get_patch(fclaw2d_patch_t* patch_2d)
{
    if(patch_2d == NULL)
    {
        return NULL;
    }
    FCLAW_ASSERT(patch_2d->user != NULL);
    return (fclaw_patch_t*) patch_2d->user;
}

void fclaw2d_patch_wrap_cb(fclaw2d_domain_t * domain_2d, 
                                 fclaw2d_patch_t * patch_2d,
                                 int blockno, int patchno, void *user)
{
    fclaw_patch_wrap_user_t* wrap = 
        (fclaw_patch_wrap_user_t*) user;
    fclaw_domain_t* domain = get_domain(domain_2d);
    fclaw_patch_t* patch = get_patch(patch_2d);
    wrap->pcb(domain, patch, blockno, patchno, wrap->user);
}

void
fclaw2d_transfer_wrap_cb(fclaw2d_domain_t * old_domain_2d,
                               fclaw2d_patch_t * old_patch_2d,
                               fclaw2d_domain_t * new_domain_2d,
                               fclaw2d_patch_t * new_patch_2d,
                               int blockno,
                               int old_patchno, int new_patchno,
                               void *user)
{
    fclaw_transfer_wrap_user_t* wrap = 
        (fclaw_transfer_wrap_user_t*) user;

    fclaw_domain_t* old_domain = get_domain(old_domain_2d);
    fclaw_patch_t* old_patch = get_patch(old_patch_2d);
    fclaw_domain_t* new_domain = get_domain(new_domain_2d);
    fclaw_patch_t* new_patch = get_patch(new_patch_2d);

    wrap->tcb(old_domain, old_patch,
              new_domain, new_patch,
              blockno, old_patchno, new_patchno,
              wrap->user);
}

void
fclaw2d_match_wrap_cb(fclaw2d_domain_t * old_domain_2d,
                            fclaw2d_patch_t * old_patch_2d,
                            fclaw2d_domain_t * new_domain_2d,
                            fclaw2d_patch_t * new_patch_2d,
                            fclaw2d_patch_relation_t newsize_2d,
                            int blockno,
                            int old_patchno, int new_patchno,
                            void *user)
{
    fclaw_match_wrap_user_t* wrap = 
        (fclaw_match_wrap_user_t*) user;

    fclaw_domain_t* old_domain = get_domain(old_domain_2d);
    fclaw_patch_t* old_patch = get_patch(old_patch_2d);
    fclaw_domain_t* new_domain = get_domain(new_domain_2d);
    fclaw_patch_t* new_patch = get_patch(new_patch_2d);

    fclaw_patch_relation_t newsize = (fclaw_patch_relation_t) newsize_2d;

    wrap->mcb(old_domain, old_patch, new_domain, new_patch, 
              newsize, blockno, old_patchno, new_patchno,
              wrap->user);
}

int 
fclaw2d_intersect_wrap_cb (fclaw2d_domain_t * domain_2d,
                        fclaw2d_patch_t * patch_2d,
                        int blockno, int patchno,
                        void *ray, double *integral,
                        void *user)
{
    fclaw_integrate_ray_wrap_user_t* wrap = 
        (fclaw_integrate_ray_wrap_user_t*) user;

    fclaw_domain_t* domain = get_domain(domain_2d);
    fclaw_patch_t* patch;

    /* if it's a leaf the user data will be set to the dim_ind patch */
    int artificial_patch = (patch_2d->user == NULL);
    if(artificial_patch)
    {
        patch = FCLAW_ALLOC_ZERO(fclaw_patch_t, 1);
        copy_patch(patch, patch_2d);
    }
    else
    {
        patch = get_patch(patch_2d);
    }

    int retval = wrap->intersect(domain, patch, blockno, patchno, ray, integral,
                                 wrap->user);
    
    if(artificial_patch)
    {
        FCLAW_FREE(patch);
        patch_2d->user = NULL;
    }

    return retval;
}

static
fclaw_domain_t* domain_wrap_meta(fclaw2d_domain_t* domain_2d)
{
    /* Only the mpi-information
       (mpicomm, mpisize and mpirank) should be copied */
    fclaw_domain_t* domain = FCLAW_ALLOC_ZERO(fclaw_domain_t, 1);
#ifndef P4_TO_P8
    domain->d2 = domain_2d;
#else
    domain->d3 = domain_2d;
#endif
    domain->refine_dim = FCLAW2D_SPACEDIM;
    domain->mpicomm = domain_2d->mpicomm;
    domain->mpisize = domain_2d->mpisize;
    domain->mpirank = domain_2d->mpirank;
    domain->attributes = domain_2d->attributes;
    domain_2d->user = domain;

    return domain;
}

int
fclaw2d_interpolate_point_wrap_cb (fclaw2d_domain_t * domain_2d,
                                fclaw2d_patch_t * patch_2d,
                                int blockno, int patchno,
                                void *point, void *user)
{
    fclaw_interpolate_point_user_wrap_t* wrap = 
        (fclaw_interpolate_point_user_wrap_t*) user;

    fclaw_domain_t* domain;

    int domain_is_meta = fclaw2d_domain_is_meta(domain_2d);
    if(domain_is_meta)
    {
        domain = domain_wrap_meta(domain_2d);
    }
    else
    {
        domain = get_domain(domain_2d);
    }

    fclaw_patch_t* patch;

    /* if it's a leaf the user data will be set to the dim_ind patch */
    int artificial_patch = (patch_2d->user == NULL);
    if(artificial_patch)
    {
        patch = FCLAW_ALLOC_ZERO(fclaw_patch_t, 1);
        copy_patch(patch, patch_2d);
    }
    else
    {
        patch = get_patch(patch_2d);
    }

    int retval = wrap->interpolate(domain, patch, blockno, patchno, point,
                                   wrap->user);

    if(artificial_patch)
    {
        FCLAW_FREE(patch);
        patch_2d->user = NULL;
    }

    if(domain_is_meta)
    {
        FCLAW_FREE(domain);
        domain_2d->user = NULL;
    }

    return retval;
}
