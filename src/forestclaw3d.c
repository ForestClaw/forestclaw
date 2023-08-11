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

#include <fclaw2d_to_3d.h>
#include "forestclaw2d.c"

const fclaw3d_patch_flags_t fclaw3d_patch_block_face_flags[6] = {
    FCLAW3D_PATCH_ON_BLOCK_FACE_0,
    FCLAW3D_PATCH_ON_BLOCK_FACE_1,
    FCLAW3D_PATCH_ON_BLOCK_FACE_2,
    FCLAW3D_PATCH_ON_BLOCK_FACE_3,
    FCLAW3D_PATCH_ON_BLOCK_FACE_4,
    FCLAW3D_PATCH_ON_BLOCK_FACE_5
};

int
fclaw3d_domain_num_edges (const fclaw_domain_t * domain)
{
    return P8EST_EDGES;
}

int
fclaw3d_patch_edge_neighbors (fclaw_domain_t * domain,
                              int blockno, int patchno, int edgeno,
                              int *rproc, int *rblockno, int *rpatchno,
                              int *redge,
                              fclaw_patch_relation_t * neighbor_size)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    p4est_locidx_t local_num, qid;
    p4est_locidx_t edgeid, cstart, cend;
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

    local_num = block->num_patches_before + patchno;
    qid = mesh->quad_to_edge[P8EST_EDGES * local_num + edgeno];

    /* We are not yet ready for general multiblock connectivities where more
     * than four blocks meet at an edge */
    if (qid >= 0)
    {
        FCLAW_ASSERT (0 <= qid);
        if (qid >= mesh->local_num_quadrants + mesh->ghost_num_quadrants)
        {
            /* This is an inter-tree (face or edge) edge neighbor */
            edgeid =
                qid - (mesh->local_num_quadrants + mesh->ghost_num_quadrants);
            FCLAW_ASSERT (edgeid < mesh->local_num_edges);
            cstart =
                fclaw2d_array_index_locidx (mesh->edge_offset, edgeid);
            cend =
                fclaw2d_array_index_locidx (mesh->edge_offset, edgeid + 1);
            if (cstart + 1 < cend)
            {
                /* At least a five-edge, which is currently not supported */
                qid = -1;
            }
            else
            {
                FCLAW_ASSERT (cstart + 1 == cend);
                qid = fclaw2d_array_index_locidx (mesh->edge_quad, cstart);
                *redge = (int)
                    *(int8_t *) sc_array_index_int (mesh->edge_edge,
                                                    (int) cstart);
                FCLAW_ASSERT (0 <= *redge && *redge < P8EST_EDGES);
            }
        }
        else
        {
            /* for intra-tree edges we take the edge is opposite */
            *redge = edgeno ^ 3;
        }
    }

    if (qid < 0)
    {
        /* The value -1 is expected for an edge on the physical boundary */
        /* Currently we also return this for five- and more-edges */
        *neighbor_size = FCLAW_PATCH_BOUNDARY;
        *redge = -1;
    }
    else
    {
        if (qid < mesh->local_num_quadrants)
        {
            /* local quadrant may be in a different tree */
            *rproc = domain->mpirank;
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
            *rproc = mesh->ghost_to_proc[qid];
            FCLAW_ASSERT (*rproc != domain->mpirank);
            q = p4est_quadrant_array_index (&ghost->ghosts, qid);
            *rblockno = (int) q->p.piggy3.which_tree;
        }
        *rpatchno = (int) qid;
        switch (q->level - block->patches[patchno].level)
        {
        case -1:
            *neighbor_size = FCLAW_PATCH_DOUBLESIZE;
            break;
        case 0:
            *neighbor_size = FCLAW_PATCH_SAMESIZE;
            break;
        case 1:
            *neighbor_size = FCLAW_PATCH_HALFSIZE;
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* *INDENT-OFF* */
        FCLAW_ASSERT (*rproc == domain->mpirank
                      || (*rpatchno >= 0
                          && *rpatchno < mesh->ghost_num_quadrants));
        FCLAW_ASSERT (*rproc != domain->mpirank
                      || (*rblockno >= 0 && *rblockno < domain->num_blocks
                          && *rpatchno >= 0
                          && *rpatchno <
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
