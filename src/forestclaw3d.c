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
fclaw3d_domain_num_edges (const fclaw3d_domain_t * domain)
{
    return P8EST_EDGES;
}

int
fclaw3d_patch_edge_neighbors (fclaw2d_domain_t * domain,
                              int blockno, int patchno, int edgeno,
                              int rproc[2], int *rblockno, int rpatchno[2],
                              int *redge,
                              fclaw2d_patch_relation_t * neighbor_size)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    p4est_locidx_t local_num, qid, qte;
    p4est_locidx_t edge_offset_i, edgeid, cstart, cend;
    const p4est_quadrant_t *q;
    p4est_tree_t *rtree;
    fclaw2d_block_t *block;
    int i, num_neighbors;

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
    qte = mesh->quad_to_edge[P8EST_EDGES * local_num + edgeno];

    /* We are not yet ready for general multiblock connectivities where more
     * than four blocks meet at an edge */
    if (qte >= 0)
    {
        /* has neighbor. process and get neighbors */
        if (qte >= mesh->local_num_quadrants + mesh->ghost_num_quadrants)
        {
            /* This is an inter-tree (face or edge) edge neighbor 
               or half/double sized neighbor */
            edge_offset_i =
                qte - (mesh->local_num_quadrants + mesh->ghost_num_quadrants);
            FCLAW_ASSERT (edge_offset_i < mesh->local_num_edges);

            cstart = fclaw2d_array_index_locidx (mesh->edge_offset, edge_offset_i);
            cend = fclaw2d_array_index_locidx (mesh->edge_offset, edge_offset_i + 1);

            /* get value in edge_edge array */
            edgeid = *(int8_t *) sc_array_index_int (mesh->edge_edge,
                                                    (int) cstart);
            if ((edgeid >= 0 && cstart + 1 < cend) || cstart + 2 < cend)
            {
                /* At least a five-edge, which is currently not supported */
                num_neighbors = 0;
                *neighbor_size = FCLAW2D_PATCH_BOUNDARY;
                *redge = -1;
            }
            else
            {
                /* at least have one neighbor, get the first neighbor */
                rpatchno[0] = fclaw2d_array_index_locidx (mesh->edge_quad, cstart);
                /* decode */
                if(edgeid < 0)
                {
                    /* half sized neighbor */
                    num_neighbors = 2;
                    *neighbor_size = FCLAW2D_PATCH_HALFSIZE;
                    *redge = (edgeid + 24)%12;
                    /* get the second neighbor */
                    rpatchno[1] = fclaw2d_array_index_locidx (mesh->edge_quad, cstart+1);
                }
                else if (edgeid > 23)
                {
                    /* double sized neighbor */
                    num_neighbors = 1;
                    *neighbor_size = FCLAW2D_PATCH_DOUBLESIZE;
                    *redge = (edgeid - 24)%12;
                }
                else 
                {
                    /* same sized neighbor inter-tree */
                    num_neighbors = 1;
                    *neighbor_size = FCLAW2D_PATCH_SAMESIZE;
                    *redge = edgeid%12;
                }
                FCLAW_ASSERT (0 <= *redge && *redge < P8EST_EDGES);
            }
        }
        else
        {
            /* for same size intra-tree edges we take the edge is opposite */
            num_neighbors = 1;
            *neighbor_size = FCLAW2D_PATCH_SAMESIZE;
            rpatchno[0] = qte;
            *redge = edgeno ^ 3;
        }
    }
    else
    {
        /* The value -1 is expected for an edge on the physical boundary */
        /* Currently we also return this for five- and more-edges */
        num_neighbors = 0;
        *neighbor_size = FCLAW2D_PATCH_BOUNDARY;
        *redge = -1;
    }

    /* get rproc and rpatchno for each neighbor */
    for(i = 0; i < num_neighbors; i++)
    {
        qid = rpatchno[i];
        if (qid < mesh->local_num_quadrants)
        {
            /* local quadrant may be in a different tree */
            rproc[i] = domain->mpirank;
            *rblockno = (int) mesh->quad_to_tree[qid];
            rtree = p4est_tree_array_index (p4est->trees,
                                            (p4est_topidx_t) *rblockno);
            FCLAW_ASSERT (rtree->quadrants_offset <= qid);
            qid -= rtree->quadrants_offset;     /* relative to tree */
            q = p4est_quadrant_array_index (&rtree->quadrants, qid);
        }
        else
        {
            qid -= mesh->local_num_quadrants;   /* relative to ghosts */
            FCLAW_ASSERT (qid < mesh->ghost_num_quadrants);
            rproc[i] = mesh->ghost_to_proc[qid];
            FCLAW_ASSERT (rproc[i] != domain->mpirank);
            q = p4est_quadrant_array_index (&ghost->ghosts, qid);
            *rblockno = (int) q->p.piggy3.which_tree;
        }
        rpatchno[i] = (int) qid;

        /* *INDENT-OFF* */
        FCLAW_ASSERT (rproc[i] == domain->mpirank
                      || (rpatchno[i] >= 0
                          && rpatchno[i] < mesh->ghost_num_quadrants));
        FCLAW_ASSERT (rproc[i] != domain->mpirank
                      || (*rblockno >= 0 && *rblockno < domain->num_blocks
                          && rpatchno[i] >= 0
                          && rpatchno[i] <
                             domain->blocks[*rblockno].num_patches));
        /* *INDENT-ON* */
    }

    return *neighbor_size != FCLAW2D_PATCH_BOUNDARY;
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
fclaw3d_patch_transform_edge (fclaw2d_patch_t * ipatch,
                              fclaw2d_patch_t * opatch,
                              int iedge, int is_block_boundary,
                              int mx, int my, int mz,
                              int based, int *i, int *j, int *k)
{
    FCLAW_ASSERT (ipatch->level == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && ipatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);
    FCLAW_ASSERT (ipatch->zlower >= 0. && ipatch->zlower < 1.);
    FCLAW_ASSERT (opatch->zlower >= 0. && opatch->zlower < 1.);

    FCLAW_ASSERT (mx >= 1 && my >= 1 && mz >= 1);
    FCLAW_ASSERT (based == 0 || based == 1);

    double Rmxmymz[P4EST_DIM], xyzshift[3];
    int shiftinds[2];
    Rmxmymz[0] = (double) (1 << ipatch->level) * (double) mx;
    Rmxmymz[1] = (double) (1 << ipatch->level) * (double) my;
    Rmxmymz[2] = (double) (1 << ipatch->level) * (double) mz;

    xyzshift[0] = xyzshift[1] = xyzshift[2] = 0.;
    if (is_block_boundary)
    {
        /* We need to add/substract the shift due to translation of one block. */
        /* determine in which dimensions we need to shift */
        if ((iedge & 12) == 0)
        {
            /* This edge is parallel to the x-axis */
            shiftinds[0] = 1;   /* first shift is in y-dimension */
            shiftinds[1] = 2;   /* second shift is in z-dimension */
        }
        else if ((iedge & 12) == 4)
        {
            /* This edge is parallel to the y-axis */
            shiftinds[0] = 0;   /* first shift is in x-dimension */
            shiftinds[1] = 2;   /* second shift is in z-dimension */
        }
        else
        {
            /* This edge is parallel to the z-axis */
            FCLAW_ASSERT ((iedge & 12) == 8);
            shiftinds[0] = 0;   /* first shift is in x-dimension */
            shiftinds[1] = 1;   /* second shift is in y-dimension */
        }

        /* shift in previously determined dimensions based on edge number */
        if ((iedge & 1) == 0)
        {
            xyzshift[shiftinds[0]] = +1.;
        }
        else
        {
            xyzshift[shiftinds[0]] = -1.;
        }
        if ((iedge & 2) == 0)
        {
            xyzshift[shiftinds[1]] = +1.;
        }
        else
        {
            xyzshift[shiftinds[1]] = -1.;
        }
        /* verify the blocks are edge-neighbors */
        FCLAW_ASSERT ((xyzshift[0] == 0.) ||
                      (xyzshift[0] == 1.
                       && fabs (opatch->xupper - 1.) < SC_1000_EPS)
                      || (xyzshift[0] == -1.
                          && fabs (opatch->xlower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[1] == 0.)
                      || (xyzshift[1] == 1.
                          && fabs (opatch->yupper - 1.) < SC_1000_EPS)
                      || (xyzshift[1] == -1.
                          && fabs (opatch->ylower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[2] == 0.)
                      || (xyzshift[2] == 1.
                          && fabs (opatch->zupper - 1.) < SC_1000_EPS)
                      || (xyzshift[2] == -1.
                          && fabs (opatch->zlower) < SC_1000_EPS));
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    *i +=
        (int) ((ipatch->xlower - opatch->xlower + xyzshift[0]) * Rmxmymz[0]);
    *j +=
        (int) ((ipatch->ylower - opatch->ylower + xyzshift[1]) * Rmxmymz[1]);
    *k +=
        (int) ((ipatch->zlower - opatch->zlower + xyzshift[2]) * Rmxmymz[2]);
}

void
fclaw3d_patch_transform_edge2 (fclaw3d_patch_t * ipatch,
                               fclaw3d_patch_t * opatch,
                               int iedge, int is_block_boundary,
                               int mx, int my, int mz, int based,
                               int i[], int j[], int k[])
{
    FCLAW_ASSERT (ipatch->level + 1 == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && opatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);
    FCLAW_ASSERT (ipatch->zlower >= 0. && ipatch->zlower < 1.);
    FCLAW_ASSERT (opatch->zlower >= 0. && opatch->zlower < 1.);

    FCLAW_ASSERT (mx >= 1 && my >= 1 && mz >= 1);
    FCLAW_ASSERT (based == 0 || based == 1);

    int kt, kn, ks;
    int di, dj, dk;
    double Rmxmymz[P4EST_DIM], xyzshift[3];
    int shiftinds[2];
    Rmxmymz[0] = (double) (1 << opatch->level) * (double) mx;
    Rmxmymz[1] = (double) (1 << opatch->level) * (double) my;
    Rmxmymz[2] = (double) (1 << opatch->level) * (double) mz;

    xyzshift[0] = xyzshift[1] = xyzshift[2] = 0.;
    if (is_block_boundary)
    {
        /* We need to add/substract the shift due to translation of one block. */
        /* determine in which dimensions we need to shift */
        if ((iedge & 12) == 0)
        {
            /* This edge is parallel to the x-axis */
            shiftinds[0] = 1;   /* first shift is in y-dimension */
            shiftinds[1] = 2;   /* second shift is in z-dimension */
        }
        else if ((iedge & 12) == 4)
        {
            /* This edge is parallel to the y-axis */
            shiftinds[0] = 0;   /* first shift is in x-dimension */
            shiftinds[1] = 2;   /* second shift is in z-dimension */
        }
        else
        {
            /* This edge is parallel to the z-axis */
            FCLAW_ASSERT ((iedge & 12) == 8);
            shiftinds[0] = 0;   /* first shift is in x-dimension */
            shiftinds[1] = 1;   /* second shift is in y-dimension */
        }

        /* shift in previously determined dimensions based on edge number */
        if ((iedge & 1) == 0)
        {
            xyzshift[shiftinds[0]] = +1.;
        }
        else
        {
            xyzshift[shiftinds[0]] = -1.;
        }
        if ((iedge & 2) == 0)
        {
            xyzshift[shiftinds[1]] = +1.;
        }
        else
        {
            xyzshift[shiftinds[1]] = -1.;
        }
        /* verify the blocks are edge-neighbors */
        FCLAW_ASSERT ((xyzshift[0] == 0.) ||
                      (xyzshift[0] == 1.
                       && fabs (opatch->xupper - 1.) < SC_1000_EPS)
                      || (xyzshift[0] == -1.
                          && fabs (opatch->xlower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[1] == 0.)
                      || (xyzshift[1] == 1.
                          && fabs (opatch->yupper - 1.) < SC_1000_EPS)
                      || (xyzshift[1] == -1.
                          && fabs (opatch->ylower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[2] == 0.)
                      || (xyzshift[2] == 1.
                          && fabs (opatch->zupper - 1.) < SC_1000_EPS)
                      || (xyzshift[2] == -1.
                          && fabs (opatch->zlower) < SC_1000_EPS));
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    di = based +
        (int) ((ipatch->xlower - opatch->xlower + xyzshift[0]) * Rmxmymz[0] +
               2. * (*i - based));
    dj = based +
        (int) ((ipatch->ylower - opatch->ylower + xyzshift[1]) * Rmxmymz[1] +
               2. * (*j - based));
    dk = based +
        (int) ((ipatch->zlower - opatch->zlower + xyzshift[2]) * Rmxmymz[2] +
               2. * (*k - based));

    /* Without any rotation, the order of child cells is canonical */
    for (ks = 0; ks < 2; ++ks)
    {
        for (kt = 0; kt < 2; ++kt)
        {
            for (kn = 0; kn < 2; ++kn)
            {
                i[4 * ks + 2 * kt + kn] = di + kn;
                j[4 * ks + 2 * kt + kn] = dj + kt;
                k[4 * ks + 2 * kt + kn] = dk + ks;
            }
        }
    }
}
