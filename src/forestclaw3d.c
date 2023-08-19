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

#define P4_TO_P8
#include "forestclaw2d.c"

fclaw3d_domain_indirect_t *
fclaw3d_domain_indirect_begin (fclaw_domain_t * domain)
{
    //DUMMY
    return NULL;
}

void
fclaw3d_domain_indirect_end (fclaw_domain_t * domain, fclaw3d_domain_indirect_t * ind)
{
    //DUMMY
}

void
fclaw3d_domain_indirect_destroy (fclaw_domain_t * domain, fclaw3d_domain_indirect_t * ind)
{
    //DUMMY
}

const fclaw3d_patch_flags_t fclaw3d_patch_block_face_flags[6] = {
    FCLAW3D_PATCH_ON_BLOCK_FACE_0,
    FCLAW3D_PATCH_ON_BLOCK_FACE_1,
    FCLAW3D_PATCH_ON_BLOCK_FACE_2,
    FCLAW3D_PATCH_ON_BLOCK_FACE_3,
    FCLAW3D_PATCH_ON_BLOCK_FACE_4,
    FCLAW3D_PATCH_ON_BLOCK_FACE_5
};

int
fclaw_patch_edge_neighbors (fclaw_domain_t * domain,
                              int blockno, int patchno, int edgeno,
                              int rprocs[], int *rblockno, int rpatchnos[],
                              int *redge,
                              fclaw_patch_relation_t * neighbor_size)
{
    if(domain->dim != 3)
    {
        fclaw_global_essentialf("fclaw_patch_edge_neighbors is only used for 3d\n");
        exit(1);
    }

    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    const p4est_quadrant_t *q;
    p4est_tree_t *rtree;
    fclaw_block_t *block;

    FCLAW_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);

    FCLAW_ASSERT (domain->pp_owned);

    FCLAW_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    FCLAW_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    FCLAW_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    FCLAW_ASSERT (0 <= patchno && patchno < block->num_patches);
    FCLAW_ASSERT (0 <= edgeno && edgeno < P8EST_EDGES);

    p4est_locidx_t local_num = block->num_patches_before + patchno;
    p4est_locidx_t quad_to_edge = mesh->quad_to_edge[P8EST_EDGES * local_num + edgeno];

    /* We are not yet ready for general multiblock connectivities where more
     * than four blocks meet at an edge */
    int num_neighbors = 0;
    p4est_locidx_t quad_ids[2];
    if (quad_to_edge >= 0)
    {
        /* has neighbor. process and get neighbors */
        if (quad_to_edge >= mesh->local_num_quadrants + mesh->ghost_num_quadrants)
        {
            /* This is an inter-tree (face or edge) edge neighbor 
               or half/double sized neighbor */
            p4est_locidx_t edge_offset_i =
                quad_to_edge - (mesh->local_num_quadrants + mesh->ghost_num_quadrants);
            FCLAW_ASSERT (edge_offset_i < mesh->local_num_edges);

            p4est_locidx_t cstart, cend;
            cstart = fclaw2d_array_index_locidx (mesh->edge_offset, edge_offset_i);
            cend = fclaw2d_array_index_locidx (mesh->edge_offset, edge_offset_i + 1);

            /* get value in edge_edge array */
            int e = *(int8_t *) sc_array_index_int (mesh->edge_edge,
                                                    (int) cstart);
            if ((e >= 0 && cstart + 1 < cend) || cstart + 2 < cend)
            {
                /* At least a five-edge, which is currently not supported */
            }
            else
            {
                /* at least have one neighbor, get the first neighbor */
                quad_ids[0] = fclaw2d_array_index_locidx (mesh->edge_quad, cstart);
                /* decode */
                if(e < 0)
                {
                    /* half sized neighbor */
                    num_neighbors = 2;
                    *neighbor_size = FCLAW_PATCH_HALFSIZE;
                    *redge = (e + 24)%12;
                    /* get the second neighbor */
                    quad_ids[1] = fclaw2d_array_index_locidx (mesh->edge_quad, cstart+1);
                }
                else if (e > 23)
                {
                    /* double sized neighbor */
                    num_neighbors = 1;
                    *neighbor_size = FCLAW_PATCH_DOUBLESIZE;
                    *redge = (e - 24)%12;
                }
                else 
                {
                    /* same sized neighbor inter-tree */
                    num_neighbors = 1;
                    *neighbor_size = FCLAW_PATCH_SAMESIZE;
                    *redge = e%12;
                }
                FCLAW_ASSERT (0 <= *redge && *redge < P8EST_EDGES);
            }
        }
        else
        {
            /* for same size intra-tree edges we take the edge is opposite */
            num_neighbors = 1;
            *neighbor_size = FCLAW_PATCH_SAMESIZE;
            quad_ids[0] = quad_to_edge;
            *redge = edgeno ^ 3;
        }
    }
    else
    {
        /* The value -1 is expected for an edge on the physical boundary */
        /* Currently we also return this for five- and more-edges */
        *neighbor_size = FCLAW_PATCH_BOUNDARY;
        *redge = -1;
    }

    /* get rproc and rpatchno for each neighbor */
    for(int i=0; i < num_neighbors; i++)
    {
        p4est_locidx_t qid = quad_ids[i];
        if (qid < mesh->local_num_quadrants)
        {
            /* local quadrant may be in a different tree */
            rprocs[i] = domain->mpirank;
            *rblockno = (int) mesh->quad_to_tree[qid];
            rtree = p4est_tree_array_index (p4est->trees,
                                            (p4est_topidx_t) * rblockno);
            FCLAW_ASSERT (rtree->quadrants_offset <= qid);
            qid -= rtree->quadrants_offset;     /* relative to tree */
            q = p4est_quadrant_array_index (&rtree->quadrants, qid);
        }
        else
        {
            qid -= mesh->local_num_quadrants;   /* relative to ghosts */
            FCLAW_ASSERT (qid < mesh->ghost_num_quadrants);
            rprocs[i] = mesh->ghost_to_proc[qid];
            FCLAW_ASSERT (rprocs[i] != domain->mpirank);
            q = p4est_quadrant_array_index (&ghost->ghosts, qid);
            *rblockno = (int) q->p.piggy3.which_tree;
        }
        rpatchnos[i] = (int) qid;

        /* *INDENT-OFF* */
        FCLAW_ASSERT (rprocs[i] == domain->mpirank
                      || (rpatchnos[0] >= 0
                          && rpatchnos[0] < mesh->ghost_num_quadrants));
        FCLAW_ASSERT (rprocs[i] != domain->mpirank
                      || (*rblockno >= 0 && *rblockno < domain->num_blocks
                          && rpatchnos[0] >= 0
                          && rpatchnos[0] <
                             domain->blocks[*rblockno].num_patches));
        /* *INDENT-ON* */
    }

    return *neighbor_size != FCLAW_PATCH_BOUNDARY;
}

void
fclaw2d_patch_edge_swap (int *edgeno, int *redgeno)
{
    int swap;

    swap = *edgeno;
    *edgeno = *redgeno;
    *redgeno = swap;
}
void
fclaw3d_patch_transform_edge (fclaw_patch_t * ipatch,
                              fclaw_patch_t * opatch,
                              const int ftransform[],
                              int mx, int my, int mz,
                              int based, int *i, int *j, int* k
                             )
{
    //TODO actually implement this
}
