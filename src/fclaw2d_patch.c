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
#include <forestclaw2d.h>
#include <p4est_base.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_domain.h>

struct fclaw2d_patch_data
{
    fclaw2d_patch_relation_t face_neighbors[4];
    fclaw2d_patch_relation_t corner_neighbors[4];
    int corners[4];
    int block_corner_count[4];
    int on_coarsefine_interface;
    int has_finegrid_neighbors;
    int neighbors_set;

    void *user_patch; /* Start of attempt to "virtualize" the user patch. */
};

fclaw2d_patch_data_t*
fclaw2d_patch_get_user_data(fclaw2d_patch_t* patch)
{
    return (fclaw2d_patch_data_t *) patch->user;
}

void*
fclaw2d_patch_get_user_patch(fclaw2d_patch_t* patch)

{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    FCLAW_ASSERT(pdata != NULL);
    return pdata->user_patch;
}


void fclaw2d_patch_set_face_type(fclaw2d_patch_t *patch,int iface,
                                 fclaw2d_patch_relation_t face_type)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    pdata->face_neighbors[iface] = face_type;
}

void fclaw2d_patch_set_corner_type(fclaw2d_patch_t *patch,int icorner,
                                   fclaw2d_patch_relation_t corner_type)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    pdata->corner_neighbors[icorner] = corner_type;
    pdata->corners[icorner] = 1;
}

void fclaw2d_patch_set_missing_corner(fclaw2d_patch_t *patch,int icorner)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    pdata->corners[icorner] = 0;
}

fclaw2d_patch_relation_t fclaw2d_patch_get_face_type(fclaw2d_patch_t* patch,
                                                     int iface)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    FCLAW_ASSERT(pdata->neighbors_set != 0);
    FCLAW_ASSERT(0 <= iface && iface < 4);
    return pdata->face_neighbors[iface];
}

fclaw2d_patch_relation_t fclaw2d_patch_get_corner_type(fclaw2d_patch_t* patch,
                                                       int icorner)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    FCLAW_ASSERT(pdata->corners[icorner] != 0);
    FCLAW_ASSERT(pdata->neighbors_set != 0);
    return pdata->corner_neighbors[icorner];
}

int fclaw2d_patch_corner_is_missing(fclaw2d_patch_t* patch,
                                    int icorner)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    return !pdata->corners[icorner];
}

void fclaw2d_patch_neighbors_set(fclaw2d_patch_t* patch)
{
    int iface, icorner;
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    FCLAW_ASSERT(pdata->neighbors_set == 0);

    pdata->has_finegrid_neighbors = 0;
    pdata->on_coarsefine_interface = 0;
    for (iface = 0; iface < 4; iface++)
    {
        fclaw2d_patch_relation_t nt;
        nt = pdata->face_neighbors[iface];
        if (nt == FCLAW2D_PATCH_HALFSIZE || (nt == FCLAW2D_PATCH_DOUBLESIZE))
        {
            pdata->on_coarsefine_interface = 1;
            if (nt == FCLAW2D_PATCH_HALFSIZE)
            {
                pdata->has_finegrid_neighbors = 1;
            }
        }
    }

    for (icorner = 0; icorner < 4; icorner++)
    {
        fclaw2d_patch_relation_t nt;
        int has_corner = pdata->corners[icorner];
        if (has_corner)
        {
            nt = pdata->corner_neighbors[icorner];
            if ((nt == FCLAW2D_PATCH_HALFSIZE) || (nt == FCLAW2D_PATCH_DOUBLESIZE))
            {
                pdata->on_coarsefine_interface = 1;
                if (nt == FCLAW2D_PATCH_HALFSIZE)
                {
                    pdata->has_finegrid_neighbors = 1;
                }
            }
        }
    }
    pdata->neighbors_set = 1;
}

void fclaw2d_patch_neighbors_reset(fclaw2d_patch_t* patch)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    pdata->neighbors_set = 0;
}

int fclaw2d_patch_neighbor_type_set(fclaw2d_patch_t* patch)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    return pdata->neighbors_set;
}


int fclaw2d_patch_has_finegrid_neighbors(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    return pdata->has_finegrid_neighbors;
}

int fclaw2d_patch_on_coarsefine_interface(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(patch);
    return pdata->on_coarsefine_interface;
}


int
fclaw2d_patch_on_parallel_boundary (const fclaw2d_patch_t * patch)
{
    return patch->flags & FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY ? 1 : 0;
}

int* fclaw2d_patch_block_corner_count(fclaw2d_domain_t* domain,
                                      fclaw2d_patch_t* this_patch)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(this_patch);
    return pdata->block_corner_count;
}

void fclaw2d_patch_set_block_corner_count(fclaw2d_domain_t* domain,
                                          fclaw2d_patch_t* this_patch,
                                          int icorner, int block_corner_count)
{
    fclaw2d_patch_data_t *pdata = fclaw2d_patch_get_user_data(this_patch);
    pdata->block_corner_count[icorner] = block_corner_count;
}


void
fclaw2d_domain_iterate_level_mthread (fclaw2d_domain_t * domain, int level,
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


/* -------------------------------------------------------
   Relies on ClawPatch
   To get the links to ClawPatch out, we need to set up
   several virtual functions so that "patch" knows how to
   interpolate, average, etc.  Right now, only ClawPatch knows
   how to do this.
   -------------------------------------------------------- */

void fclaw2d_patch_data_new(fclaw2d_domain_t* domain,
                            fclaw2d_patch_t* this_patch)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);

    /* Initialize user data */
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    this_patch->user = (void *) pdata;

    /* create new user data */
    pdata->user_patch = patch_vt.patch_new();
    ++ddata->count_set_clawpatch;
    pdata->neighbors_set = 0;
}

void fclaw2d_patch_data_delete(fclaw2d_domain_t* domain,
                               fclaw2d_patch_t *this_patch)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t*) this_patch->user;

    if (pdata != NULL)
    {
        fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
        patch_vt.patch_delete(pdata->user_patch);
        ++ddata->count_delete_clawpatch;

        FCLAW2D_FREE(pdata);
        this_patch->user = NULL;
    }
}

void fclaw2d_set_patch_vtable(fclaw2d_domain_t* domain, fclaw2d_patch_vtable_t *patch_vt)
{
    fclaw2d_domain_attribute_add (domain,"patch_vtable",patch_vt);
}

fclaw2d_patch_vtable_t fclaw2d_get_patch_vtable(fclaw2d_domain_t* domain)
{
    fclaw2d_patch_vtable_t *patch_vt;
    patch_vt = (fclaw2d_patch_vtable_t*) fclaw2d_domain_attribute_access(domain,"patch_vtable",NULL);
    FCLAW_ASSERT(patch_vt != NULL);
    return *patch_vt;
}

void fclaw2d_patch_pack_local_ghost(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              double *patch_data,
                              int time_interp)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.ghost_pack(domain, 
                        this_patch,
                        patch_data,
                        time_interp);
}

void fclaw2d_patch_build_remote_ghost(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int blockno,
                               int patchno,
                               void *user)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.build_ghost(domain,this_patch,blockno,
                         patchno,(void*) user);
}

void fclaw2d_patch_unpack_remote_ghost(fclaw2d_domain_t* domain,
                                       fclaw2d_patch_t* this_patch,
                                       int this_block_idx,
                                       int this_patch_idx,
                                       double *qdata, fclaw_bool time_interp)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.ghost_unpack(domain, this_patch, this_block_idx, 
                          this_patch_idx, qdata, time_interp);
}

size_t fclaw2d_patch_ghost_packsize(fclaw2d_domain_t* domain)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    return patch_vt.ghost_packsize(domain);
}

void fclaw2d_patch_alloc_local_ghost(fclaw2d_domain_t* domain,
                                     fclaw2d_patch_t* this_patch,
                                     void** q)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.local_ghost_alloc(domain, this_patch, q);
}

void fclaw2d_patch_free_local_ghost(fclaw2d_domain_t* domain,
                                    void **q)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.local_ghost_free(domain, q);
}

void cb_fclaw2d_patch_partition_pack(fclaw2d_domain_t *domain,
                                     fclaw2d_patch_t *this_patch,
                                     int this_block_idx,
                                     int this_patch_idx,
                                     void *user)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.partition_pack(domain,
                            this_patch,
                            this_block_idx,
                            this_patch_idx,
                            user);
}

void cb_fclaw2d_patch_partition_unpack(fclaw2d_domain_t *domain,
                                     fclaw2d_patch_t *this_patch,
                                     int this_block_idx,
                                     int this_patch_idx,
                                     void *user)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    patch_vt.partition_unpack(domain,
                            this_patch,
                            this_block_idx,
                            this_patch_idx,
                            user);
}

size_t fclaw2d_patch_partition_packsize(fclaw2d_domain_t* domain)
{
    fclaw2d_patch_vtable_t patch_vt = fclaw2d_get_patch_vtable(domain);
    return patch_vt.partition_packsize(domain);
}