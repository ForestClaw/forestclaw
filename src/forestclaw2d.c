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

#include "forestclaw2d.h"
#include <p4est_bits.h>

int
fclaw2d_domain_dimension (fclaw2d_domain_t * domain)
{
    return P4EST_DIM;           /* space dimension */
}

int
fclaw2d_domain_num_faces (fclaw2d_domain_t * domain)
{
    return P4EST_FACES;         /* 2 * DIM; number of cube faces */
}

int
fclaw2d_domain_num_corners (fclaw2d_domain_t * domain)
{
    return P4EST_CHILDREN;      /* 2 ** DIM; number of cube corners */
}

int
fclaw2d_domain_num_face_corners (fclaw2d_domain_t * domain)
{
    return P4EST_HALF;          /* 2 ** (DIM - 1); corners per face */
}

int
fclaw2d_domain_num_orientations (fclaw2d_domain_t * domain)
{
    return P4EST_FACES * P4EST_HALF;
}

void
fclaw2d_domain_corner_faces (fclaw2d_domain_t * domain,
                             int icorner, int faces[2])
{
    P4EST_ASSERT (0 <= icorner && icorner < P4EST_CHILDREN);
    faces[0] = p4est_corner_faces[icorner][0];
    faces[1] = p4est_corner_faces[icorner][1];
}

void
fclaw2d_domain_iterate_level (fclaw2d_domain_t * domain, int level,
                              fclaw2d_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    P4EST_ASSERT (0 <= level && level <= domain->possible_maxlevel);
    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        for (patch = block->patchbylevel[level];
             patch != NULL; patch = patch->next)
        {
            j = (int) (patch - block->patches);
            P4EST_ASSERT (0 <= j && j < block->num_patches);
            P4EST_ASSERT (patch->level == level);
            pcb (domain, patch, i, j, user);
        }
    }
}

void
fclaw2d_domain_iterate_patches (fclaw2d_domain_t * domain,
                                fclaw2d_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    for (i = 0; i < domain->num_blocks; i++)
    {
        block = domain->blocks + i;
        for (j = 0; j < block->num_patches; j++)
        {
            patch = block->patches + j;
            pcb (domain, patch, i, j, user);
        }
    }
}

int
fclaw2d_patch_boundary_type (fclaw2d_domain_t * domain,
                             int blockno, int patchno,
                             int boundaries[P4EST_FACES])
{
    int faceno;
    int anyboundary;
    int8_t qtf;
    p4est_locidx_t totalleaf;
    p4est_locidx_t qtq;
    p4est_t *p4est;
    p4est_tree_t *tree;
    p4est_mesh_t *mesh;
    fclaw2d_block_t *block;

    P4EST_ASSERT (domain->pp_owned);
    anyboundary = 0;
    p4est = domain->pp->p4est;
    mesh = domain->pp->match_aux ? domain->pp->mesh_aux : domain->pp->mesh;

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    totalleaf = tree->quadrants_offset + (p4est_locidx_t) patchno;
    P4EST_ASSERT (0 <= totalleaf && totalleaf < p4est->local_num_quadrants);
    for (faceno = 0; faceno < P4EST_FACES; ++faceno)
    {
        qtq = mesh->quad_to_quad[P4EST_FACES * totalleaf + faceno];
        qtf = mesh->quad_to_face[P4EST_FACES * totalleaf + faceno];
        if (qtq == totalleaf && qtf == faceno)
        {
            P4EST_ASSERT (block->is_boundary[faceno]);
            anyboundary = 1;
            boundaries[faceno] = 1;
        }
        else
        {
            boundaries[faceno] = 0;
        }
    }

    return anyboundary;
}

static void
fclaw2d_patch_encode_neighbor (fclaw2d_domain_t * domain, p4est_mesh_t * mesh,
                               p4est_locidx_t qtq, int *proc, int *blockno,
                               int *patchno)
{
    p4est_quadrant_t *ghost;

    P4EST_ASSERT (domain->pp_owned);
    P4EST_ASSERT (0 <= qtq);
    P4EST_ASSERT (qtq <
                  mesh->local_num_quadrants + mesh->ghost_num_quadrants);
    if (qtq < mesh->local_num_quadrants)
    {
        /* processor-local neighbor */
        *proc = domain->mpirank;
        *blockno = domain->patch_to_block[qtq];
        *patchno = (int) qtq;   /* patch number within the block as usual */
    }
    else
    {
        /* off-processor ghost neighbor */
        qtq -= mesh->local_num_quadrants;
        *proc = mesh->ghost_to_proc[qtq];
        ghost = p4est_quadrant_array_index (&domain->pp->ghost->ghosts, qtq);
        P4EST_ASSERT (0 <= ghost->p.piggy3.which_tree);
        P4EST_ASSERT (ghost->p.piggy3.which_tree <
                      domain->pp->conn->num_trees);
        *blockno = (int) ghost->p.piggy3.which_tree;
        *patchno = (int) ghost->p.piggy3.local_num;
        P4EST_ASSERT (*patchno == (int) mesh->ghost_to_index[qtq]);
    }
}

fclaw2d_face_neighbor_t
fclaw2d_patch_face_neighbors (fclaw2d_domain_t * domain,
                              int blockno, int patchno, int faceno,
                              int rproc[P4EST_HALF], int *rblockno,
                              int rpatchno[P4EST_HALF], int *rfaceno)
{
    const int num_orientations = fclaw2d_domain_num_orientations (domain);
    int k;
    int hblockno[P4EST_HALF];
    int8_t qtf;
    p4est_locidx_t totalleaf;
    p4est_locidx_t qtq, *qth;
    p4est_t *p4est;
    p4est_tree_t *tree;
    p4est_mesh_t *mesh;
    fclaw2d_block_t *block;

    P4EST_ASSERT (domain->pp_owned);
    p4est = domain->pp->p4est;
    mesh = domain->pp->match_aux ? domain->pp->mesh_aux : domain->pp->mesh;

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    totalleaf = tree->quadrants_offset + (p4est_locidx_t) patchno;
    P4EST_ASSERT (0 <= totalleaf && totalleaf < p4est->local_num_quadrants);
    P4EST_ASSERT (0 <= faceno && faceno < P4EST_FACES);

    qtq = mesh->quad_to_quad[P4EST_FACES * totalleaf + faceno];
    qtf = mesh->quad_to_face[P4EST_FACES * totalleaf + faceno];
    if (qtq == totalleaf && qtf == faceno)
    {
        /* physical domain boundary encoded by same patch face */
        P4EST_ASSERT (block->is_boundary[faceno]);
        rproc[0] = domain->mpirank;
        *rblockno = blockno;
        rpatchno[0] = patchno;
        *rfaceno = faceno;
        for (k = 1; k < P4EST_HALF; ++k)
        {
            rproc[k] = -1;
            rpatchno[k] = -1;
        }
        return FCLAW2D_FACE_NEIGHBOR_BOUNDARY;
    }
    else if (qtf < 0)
    {
        /* half-size face neighbors */
        qth = (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, qtq);
        for (k = 0; k < P4EST_HALF; ++k)
        {
            fclaw2d_patch_encode_neighbor (domain, mesh, qth[k],
                                           rproc + k, hblockno + k,
                                           rpatchno + k);
            P4EST_ASSERT (k == 0 || hblockno[k - 1] == hblockno[k]);
        }
        *rblockno = hblockno[0];
        *rfaceno = qtf + num_orientations;
        P4EST_ASSERT (*rfaceno >= 0);
        return FCLAW2D_FACE_NEIGHBOR_HALFSIZE;
    }
    else
    {
        /* one same-size or double-size neighbor */
        fclaw2d_patch_encode_neighbor (domain, mesh, qtq,
                                       rproc, rblockno + 0, rpatchno);
        for (k = 1; k < P4EST_HALF; ++k)
        {
            rproc[k] = -1;
            rpatchno[k] = -1;
        }
        if (qtf < num_orientations)
        {
            /* same-size neighbor */
            *rfaceno = (int) qtf;
            P4EST_ASSERT (0 <= *rfaceno && *rfaceno < num_orientations);
            return FCLAW2D_FACE_NEIGHBOR_SAMESIZE;
        }
        else
        {
            /* double-size neighbor */
            *rfaceno = (int) qtf % num_orientations;
            /* the number of our patch within the bigger neighbor subfaces */
            rproc[1] = (int) qtf / num_orientations - 1;
            P4EST_ASSERT (0 <= rproc[1] && rproc[1] < P4EST_HALF);
            return FCLAW2D_FACE_NEIGHBOR_DOUBLESIZE;
        }
    }
}

int
fclaw2d_patch_corner_neighbors (fclaw2d_domain_t * domain,
                                int blockno, int patchno, int cornerno)
{
    const p4est_quadrant_t *q;
    p4est_quadrant_t r;
    p4est_t *p4est;
    p4est_tree_t *tree;
    fclaw2d_block_t *block;

    P4EST_ASSERT (domain->pp_owned);
    p4est = domain->pp->p4est;

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    q = p4est_quadrant_array_index (&tree->quadrants, patchno);

    P4EST_ASSERT (0 <= cornerno && cornerno < P4EST_CHILDREN);
    p4est_quadrant_corner_neighbor (q, cornerno, &r);

    return
        (int) sc_array_bsearch (&tree->quadrants, &r, p4est_quadrant_compare);
}

void
fclaw2d_patch_mark_refine (fclaw2d_domain_t * domain, int blockno,
                           int patchno)
{
    int totalpatchno;
    fclaw2d_block_t *block;

    P4EST_ASSERT (domain->pp_owned);
    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);
    totalpatchno = block->num_patches_before + patchno;
    P4EST_ASSERT (0 <= totalpatchno
                  && totalpatchno < domain->num_patches_all);

    domain->pp->flags[totalpatchno] = P4EST_WRAP_REFINE;
}

void
fclaw2d_patch_mark_coarsen (fclaw2d_domain_t * domain, int blockno,
                            int patchno)
{
    int totalpatchno;
    fclaw2d_block_t *block;

    P4EST_ASSERT (domain->pp_owned);
    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);
    totalpatchno = block->num_patches_before + patchno;
    P4EST_ASSERT (0 <= totalpatchno
                  && totalpatchno < domain->num_patches_all);

    domain->pp->flags[totalpatchno] = P4EST_WRAP_COARSEN;
}

/* Iterate over patches at level 'level' that didn't change upon regridding */
/* 'level' here refers to the level of the old patch */
void
fclaw2d_domain_iterate_unchanged (fclaw2d_domain_t * old_domain,
                                  fclaw2d_domain_t * new_domain, int level,
                                  fclaw2d_match_unchanged_callback_t cb_user,
                                  void *user)
{
}

/* Iterate over patches which have been refined */
void
fclaw2d_domain_iterate_refined (fclaw2d_domain_t * old_domain,
                                fclaw2d_domain_t * new_domain, int level,
                                fclaw2d_match_refined_callback_t cb_user,
                                void *user)
{
}

/* Iterate over patches which have been coarsened */
void
fclaw2d_domain_iterate_coarsened (fclaw2d_domain_t * old_domain,
                                  fclaw2d_domain_t * new_domain, int level,
                                  fclaw2d_match_coarsened_callback_t cb_user,
                                  void *user)
{
}
