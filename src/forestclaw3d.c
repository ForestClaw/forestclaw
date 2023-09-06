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
fclaw_patch_edge_neighbors (fclaw_domain_t * domain,
                              int blockno, int patchno, int edgeno,
                              int rprocs_out[], int *rblockno_out,
                              int rpatchnos_out[],
                              int *redge_out,
                              fclaw_patch_relation_t * neighbor_size_out)
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
    fclaw_patch_relation_t neighbor_size;
    int redge;
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
                    neighbor_size = FCLAW_PATCH_HALFSIZE;
                    redge = (e + 24)%12;
                    /* get the second neighbor */
                    quad_ids[1] = fclaw2d_array_index_locidx (mesh->edge_quad, cstart+1);
                }
                else if (e > 23)
                {
                    /* double sized neighbor */
                    num_neighbors = 1;
                    neighbor_size = FCLAW_PATCH_DOUBLESIZE;
                    redge = (e - 24)%12;
                }
                else 
                {
                    /* same sized neighbor inter-tree */
                    num_neighbors = 1;
                    neighbor_size = FCLAW_PATCH_SAMESIZE;
                    redge = e%12;
                }
                FCLAW_ASSERT (0 <= redge && redge < P8EST_EDGES);
            }
        }
        else
        {
            /* for same size intra-tree edges we take the edge is opposite */
            num_neighbors = 1;
            neighbor_size = FCLAW_PATCH_SAMESIZE;
            quad_ids[0] = quad_to_edge;
            redge = edgeno ^ 3;
        }
    }
    else
    {
        /* The value -1 is expected for an edge on the physical boundary */
        /* Currently we also return this for five- and more-edges */
        neighbor_size = FCLAW_PATCH_BOUNDARY;
        redge = -1;
    }

    /* get rproc and rpatchno for each neighbor */
    int rprocs[2];
    int rblockno;
    int rpatchnos[2];
    for(int i=0; i < num_neighbors; i++)
    {
        p4est_locidx_t qid = quad_ids[i];
        if (qid < mesh->local_num_quadrants)
        {
            /* local quadrant may be in a different tree */
            rprocs[i] = domain->mpirank;
            rblockno = (int) mesh->quad_to_tree[qid];
            rtree = p4est_tree_array_index (p4est->trees,
                                            (p4est_topidx_t) rblockno);
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
            rblockno = (int) q->p.piggy3.which_tree;
        }
        rpatchnos[i] = (int) qid;

        /* *INDENT-OFF* */
        FCLAW_ASSERT (rprocs[i] == domain->mpirank
                      || (rpatchnos[0] >= 0
                          && rpatchnos[0] < mesh->ghost_num_quadrants));
        FCLAW_ASSERT (rprocs[i] != domain->mpirank
                      || (rblockno >= 0 && rblockno < domain->num_blocks
                          && rpatchnos[0] >= 0
                          && rpatchnos[0] <
                             domain->blocks[rblockno].num_patches));
        /* *INDENT-ON* */
    }

    /* set output variables
       compiler may warn about uninitialized variables, 
       and Valgrind will warn about uninitialized access.
       makes debugging easier */
    *rblockno_out = rblockno;
    *redge_out = redge;
    *neighbor_size_out = neighbor_size;
    for(int i=0; i < num_neighbors; i++)
    {
        rprocs_out[i] = rprocs[i];
        rpatchnos_out[i] = rpatchnos[i];
    }

    return neighbor_size != FCLAW_PATCH_BOUNDARY;
}

void
fclaw_patch_edge_swap (int *edgeno, int *redgeno)
{
    int swap;

    swap = *edgeno;
    *edgeno = *redgeno;
    *redgeno = swap;
}

void
fclaw3d_patch_transform_edge (fclaw_patch_t * ipatch,
                              fclaw_patch_t * opatch,
                              int iedge, int is_block_boundary,
                              int mx, int my, int mz,
                              int based, int *i, int *j, int *k)
{
    FCLAW_ASSERT (ipatch->dim == 3);
    FCLAW_ASSERT (opatch->dim == 3);
    FCLAW_ASSERT (ipatch->level == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && ipatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->d3->xlower >= 0. && ipatch->d3->xlower < 1.);
    FCLAW_ASSERT (opatch->d3->xlower >= 0. && opatch->d3->xlower < 1.);
    FCLAW_ASSERT (ipatch->d3->ylower >= 0. && ipatch->d3->ylower < 1.);
    FCLAW_ASSERT (opatch->d3->ylower >= 0. && opatch->d3->ylower < 1.);
    FCLAW_ASSERT (ipatch->d3->zlower >= 0. && ipatch->d3->zlower < 1.);
    FCLAW_ASSERT (opatch->d3->zlower >= 0. && opatch->d3->zlower < 1.);

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
                       && fabs (opatch->d3->xupper - 1.) < SC_1000_EPS)
                      || (xyzshift[0] == -1.
                          && fabs (opatch->d3->xlower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[1] == 0.)
                      || (xyzshift[1] == 1.
                          && fabs (opatch->d3->yupper - 1.) < SC_1000_EPS)
                      || (xyzshift[1] == -1.
                          && fabs (opatch->d3->ylower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[2] == 0.)
                      || (xyzshift[2] == 1.
                          && fabs (opatch->d3->zupper - 1.) < SC_1000_EPS)
                      || (xyzshift[2] == -1.
                          && fabs (opatch->d3->zlower) < SC_1000_EPS));
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    *i +=
        (int) ((ipatch->d3->xlower - opatch->d3->xlower + xyzshift[0]) * Rmxmymz[0]);
    *j +=
        (int) ((ipatch->d3->ylower - opatch->d3->ylower + xyzshift[1]) * Rmxmymz[1]);
    *k +=
        (int) ((ipatch->d3->zlower - opatch->d3->zlower + xyzshift[2]) * Rmxmymz[2]);
}

void
fclaw3d_patch_transform_edge2 (fclaw_patch_t * ipatch,
                               fclaw_patch_t * opatch,
                               int iedge, int is_block_boundary,
                               int mx, int my, int mz, int based,
                               int i[], int j[], int k[])
{
    FCLAW_ASSERT (ipatch->dim == 3);
    FCLAW_ASSERT (opatch->dim == 3);
    FCLAW_ASSERT (ipatch->level + 1 == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && opatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->d3->xlower >= 0. && ipatch->d3->xlower < 1.);
    FCLAW_ASSERT (opatch->d3->xlower >= 0. && opatch->d3->xlower < 1.);
    FCLAW_ASSERT (ipatch->d3->ylower >= 0. && ipatch->d3->ylower < 1.);
    FCLAW_ASSERT (opatch->d3->ylower >= 0. && opatch->d3->ylower < 1.);
    FCLAW_ASSERT (ipatch->d3->zlower >= 0. && ipatch->d3->zlower < 1.);
    FCLAW_ASSERT (opatch->d3->zlower >= 0. && opatch->d3->zlower < 1.);

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
                       && fabs (opatch->d3->xupper - 1.) < SC_1000_EPS)
                      || (xyzshift[0] == -1.
                          && fabs (opatch->d3->xlower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[1] == 0.)
                      || (xyzshift[1] == 1.
                          && fabs (opatch->d3->yupper - 1.) < SC_1000_EPS)
                      || (xyzshift[1] == -1.
                          && fabs (opatch->d3->ylower) < SC_1000_EPS));
        FCLAW_ASSERT ((xyzshift[2] == 0.)
                      || (xyzshift[2] == 1.
                          && fabs (opatch->d3->zupper - 1.) < SC_1000_EPS)
                      || (xyzshift[2] == -1.
                          && fabs (opatch->d3->zlower) < SC_1000_EPS));
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    di = based +
        (int) ((ipatch->d3->xlower - opatch->d3->xlower + xyzshift[0]) * Rmxmymz[0] +
               2. * (*i - based));
    dj = based +
        (int) ((ipatch->d3->ylower - opatch->d3->ylower + xyzshift[1]) * Rmxmymz[1] +
               2. * (*j - based));
    dk = based +
        (int) ((ipatch->d3->zlower - opatch->d3->zlower + xyzshift[2]) * Rmxmymz[2] +
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
