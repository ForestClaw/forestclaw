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
#include <p4est_wrap.h>

#define FCLAW2D_DOMAIN_TAG_SERIALIZE 4526

double
fclaw2d_domain_global_maximum (fclaw2d_domain_t * domain, double d)
{
    int mpiret;
    double gd;

    mpiret = MPI_Allreduce (&d, &gd, 1, MPI_DOUBLE, MPI_MAX, domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return gd;
}

double
fclaw2d_domain_global_sum (fclaw2d_domain_t * domain, double d)
{
    int mpiret;
    double gd;

    mpiret = MPI_Allreduce (&d, &gd, 1, MPI_DOUBLE, MPI_SUM, domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return gd;
}

void
fclaw2d_domain_barrier (fclaw2d_domain_t * domain)
{
    int mpiret;

    mpiret = MPI_Barrier (domain->mpicomm);
    SC_CHECK_MPI (mpiret);
}

int
fclaw2d_domain_dimension (const fclaw2d_domain_t * domain)
{
    return P4EST_DIM;           /* space dimension */
}

int
fclaw2d_domain_num_faces (const fclaw2d_domain_t * domain)
{
    return P4EST_FACES;         /* 2 * DIM; number of cube faces */
}

int
fclaw2d_domain_num_corners (const fclaw2d_domain_t * domain)
{
    return P4EST_CHILDREN;      /* 2 ** DIM; number of cube corners */
}

int
fclaw2d_domain_num_face_corners (const fclaw2d_domain_t * domain)
{
    return P4EST_HALF;          /* 2 ** (DIM - 1); corners per face */
}

int
fclaw2d_domain_num_orientations (const fclaw2d_domain_t * domain)
{
    return P4EST_FACES * P4EST_HALF;
}

void
fclaw2d_domain_corner_faces (const fclaw2d_domain_t * domain,
                             int icorner, int faces[2])
{
    P4EST_ASSERT (0 <= icorner && icorner < P4EST_CHILDREN);
    faces[0] = p4est_corner_faces[icorner][0];
    faces[1] = p4est_corner_faces[icorner][1];
}

int
fclaw2d_patch_corner_dimension (const fclaw2d_patch_t * patch, int cornerno)
{
    const int childid = fclaw2d_patch_childid (patch);

    P4EST_ASSERT (0 <= cornerno && cornerno < P4EST_CHILDREN);

    return (patch->level == 0 ||
            cornerno == childid ||
            cornerno == P4EST_CHILDREN - 1 - childid) ? 0 : 1;
}

int
fclaw2d_patch_childid (const fclaw2d_patch_t * patch)
{
    const int childid = patch->flags & FCLAW2D_PATCH_CHILDID;

    P4EST_ASSERT (0 <= childid && childid < P4EST_CHILDREN);

    return childid;
}

int
fclaw2d_patch_is_first_sibling (const fclaw2d_patch_t * patch)
{
    return patch->flags & FCLAW2D_PATCH_FIRST_SIBLING ? 1 : 0;
}

int
fclaw2d_patch_is_ghost (const fclaw2d_patch_t * patch)
{
    return patch->flags & FCLAW2D_PATCH_IS_GHOST ? 1 : 0;
}

void *
fclaw2d_alloc (size_t size)
{
    return sc_malloc (p4est_package_id, size);
}

void *
fclaw2d_calloc (size_t nmemb, size_t size)
{
    return sc_calloc (p4est_package_id, nmemb, size);
}

void *
fclaw2d_realloc (void *ptr, size_t size)
{
    return sc_realloc (p4est_package_id, ptr, size);
}

void
fclaw2d_free (void *ptr)
{
    sc_free (p4est_package_id, ptr);
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
             patch != NULL; patch = patch->u.next)
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

void
fclaw2d_domain_iterate_families (fclaw2d_domain_t * domain,
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
            if (fclaw2d_patch_is_first_sibling (patch))
            {
#ifdef P4EST_DEBUG
                int k;
                for (k = 0; k < P4EST_CHILDREN; ++k)
                {
                    P4EST_ASSERT (j + k < block->num_patches);
                    P4EST_ASSERT (fclaw2d_patch_childid (patch + k) == k);
                }
#endif
                pcb (domain, patch, i, j, user);
                j += P4EST_CHILDREN - 1;
            }
        }
    }
}

int
fclaw2d_patch_boundary_type (fclaw2d_domain_t * domain,
                             int blockno, int patchno,
                             int boundaries[P4EST_FACES])
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    int faceno;
    int anyboundary;
    int8_t qtf;
    p4est_locidx_t totalleaf;
    p4est_locidx_t qtq;
    p4est_tree_t *tree;
#ifdef P4EST_DEBUG
    fclaw2d_block_t *block;
#endif

    P4EST_ASSERT (domain->pp_owned);
    anyboundary = 0;

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

#ifdef P4EST_DEBUG
    block = domain->blocks + blockno;
#endif
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
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
    p4est_quadrant_t *gq;
    fclaw2d_block_t *block;

    P4EST_ASSERT (domain->pp_owned);
    P4EST_ASSERT (0 <= qtq);
    P4EST_ASSERT (qtq <
                  mesh->local_num_quadrants + mesh->ghost_num_quadrants);
    if (qtq < mesh->local_num_quadrants)
    {
        /* processor-local neighbor */
        *proc = domain->mpirank;
        *blockno = (int) mesh->quad_to_tree[qtq];
        P4EST_ASSERT ((int) wrap->p4est->first_local_tree <= *blockno);
        P4EST_ASSERT (*blockno <= (int) wrap->p4est->last_local_tree);
        block = domain->blocks + *blockno;
        qtq -= block->num_patches_before;
        P4EST_ASSERT (0 <= qtq && qtq < block->num_patches);
        *patchno = (int) qtq;   /* patch number within the block as usual */
    }
    else
    {
        /* off-processor ghost neighbor */
        qtq -= mesh->local_num_quadrants;
        P4EST_ASSERT (qtq >= 0 && qtq < mesh->ghost_num_quadrants);
        *proc = mesh->ghost_to_proc[qtq];
        gq = p4est_quadrant_array_index (&ghost->ghosts, qtq);
        P4EST_ASSERT (0 <= gq->p.piggy3.which_tree);
        P4EST_ASSERT (gq->p.piggy3.which_tree < wrap->conn->num_trees);
        *blockno = (int) gq->p.piggy3.which_tree;
        *patchno = (int) qtq;
    }
}

fclaw2d_patch_relation_t
fclaw2d_patch_face_neighbors (fclaw2d_domain_t * domain,
                              int blockno, int patchno, int faceno,
                              int rproc[P4EST_HALF], int *rblockno,
                              int rpatchno[P4EST_HALF], int *rfaceno)
{
    const int num_orientations = fclaw2d_domain_num_orientations (domain);
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    int k;
    int hblockno[P4EST_HALF];
    int8_t qtf;
    p4est_locidx_t totalleaf;
    p4est_locidx_t qtq, *qth;
    p4est_tree_t *tree;
#ifdef P4EST_DEBUG
    fclaw2d_block_t *block;
#endif

    P4EST_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);

    P4EST_ASSERT (domain->pp_owned);

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

#ifdef P4EST_DEBUG
    block = domain->blocks + blockno;
#endif
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
        return FCLAW2D_PATCH_BOUNDARY;
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
        return FCLAW2D_PATCH_HALFSIZE;
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
            return FCLAW2D_PATCH_SAMESIZE;
        }
        else
        {
            /* double-size neighbor */
            *rfaceno = (int) qtf % num_orientations;
            /* the number of our patch within the bigger neighbor subfaces */
            rproc[1] = (int) qtf / num_orientations - 1;
            P4EST_ASSERT (0 <= rproc[1] && rproc[1] < P4EST_HALF);
            return FCLAW2D_PATCH_DOUBLESIZE;
        }
    }
}

void
fclaw2d_patch_face_transformation (int faceno, int rfaceno, int ftransform[])
{
    p4est_expand_face_transform (faceno, rfaceno, ftransform);
}

#if 0

static int
fclaw2d_patch_corner_neighbors_old (fclaw2d_domain_t * domain,
                                    int blockno, int patchno, int cornerno,
                                    int *rproc, int *rblockno, int *rpatchno,
                                    fclaw2d_patch_relation_t * neighbor_size)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
#ifdef P4EST_DEBUG
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
#endif
    const p4est_topidx_t nt = (p4est_topidx_t) blockno;
    const p4est_quadrant_t *q;
    p4est_quadrant_t r, *rq;
    p4est_tree_t *tree;
    sc_array_t searr, *earr = &searr;
    sc_array_t srparr, *rparr = &srparr;
    sc_array_t sqarr, *qarr = &sqarr;
    fclaw2d_block_t *block;
    fclaw2d_patch_relation_t prel;

    P4EST_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);

    P4EST_ASSERT (domain->pp_owned);

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    q = p4est_quadrant_array_index (&tree->quadrants, patchno);
    sc_array_init (earr, sizeof (int));
    sc_array_init (rparr, sizeof (int));
    sc_array_init (qarr, sizeof (p4est_quadrant_t));

    /* TODO: Extend this to the multi-block case.
     * Then neighbors on multiple levels may exist simultaneously.
     * The if-else construct below must be reworked. */

    if (p4est_quadrant_corner_neighbor (q, cornerno, &r),
        p4est_quadrant_exists (p4est, ghost, nt, &r, earr, rparr, qarr))
    {
        prel = FCLAW2D_PATCH_SAMESIZE;
    }
    else if (q->level < P4EST_QMAXLEVEL &&
             (p4est_quadrant_half_corner_neighbor (q, cornerno, &r),
              p4est_quadrant_exists (p4est, ghost, nt, &r, earr, rparr,
                                     qarr)))
    {
        prel = FCLAW2D_PATCH_HALFSIZE;
    }
    else if (q->level > 0 && p4est_quadrant_child_id (q) == cornerno &&
             (p4est_quadrant_parent (q, &r),
              p4est_quadrant_corner_neighbor (&r, cornerno, &r),
              p4est_quadrant_exists (p4est, ghost, nt, &r, earr, rparr,
                                     qarr)))
    {
        prel = FCLAW2D_PATCH_DOUBLESIZE;
    }
    else
    {
        prel = FCLAW2D_PATCH_BOUNDARY;
    }
    *neighbor_size = prel;

    if (prel != FCLAW2D_PATCH_BOUNDARY)
    {
        P4EST_ASSERT (rparr->elem_count == 1);
        P4EST_ASSERT (qarr->elem_count == 1);
        *rproc = *(int *) sc_array_index (rparr, 0);
        rq = p4est_quadrant_array_index (qarr, 0);
        *rblockno = (p4est_topidx_t) rq->p.piggy3.which_tree;
        *rpatchno = (p4est_topidx_t) rq->p.piggy3.local_num;    /* ghost index */
        P4EST_ASSERT (*rproc == domain->mpirank ||
                      (*rpatchno >= 0
                       && *rpatchno < mesh->ghost_num_quadrants));
        P4EST_ASSERT (*rproc != domain->mpirank
                      || (*rblockno >= 0 && *rblockno < domain->num_blocks
                          && *rpatchno >= 0 && *rpatchno <
                          domain->blocks[*rblockno].num_patches));
    }

    sc_array_reset (earr);
    sc_array_reset (rparr);
    sc_array_reset (qarr);

    return prel != FCLAW2D_PATCH_BOUNDARY;
}

#endif

int
fclaw2d_patch_corner_neighbors (fclaw2d_domain_t * domain,
                                int blockno, int patchno, int cornerno,
                                int *rproc, int *rblockno, int *rpatchno,
                                fclaw2d_patch_relation_t * neighbor_size)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
#ifdef P4EST_DEBUG
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
#endif
    const p4est_topidx_t nt = (p4est_topidx_t) blockno;
    const p4est_quadrant_t *q;
    p4est_quadrant_t r, *rq;
    p4est_tree_t *tree;
    sc_array_t searr, *earr = &searr;
    sc_array_t srparr, *rparr = &srparr;
    sc_array_t sqarr, *qarr = &sqarr;
    fclaw2d_block_t *block;
    fclaw2d_patch_relation_t prel;

    P4EST_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);

    P4EST_ASSERT (domain->pp_owned);

    P4EST_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    P4EST_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    P4EST_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    P4EST_ASSERT (0 <= patchno && patchno < block->num_patches);

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    q = p4est_quadrant_array_index (&tree->quadrants, patchno);
    sc_array_init (earr, sizeof (int));
    sc_array_init (rparr, sizeof (int));
    sc_array_init (qarr, sizeof (p4est_quadrant_t));

    /* TODO: Extend this to the multi-block case.
     * Then neighbors on multiple levels may exist simultaneously.
     * The if-else construct below must be reworked. */

    if (p4est_quadrant_corner_neighbor (q, cornerno, &r),
        p4est_quadrant_exists (p4est, ghost, nt, &r, earr, rparr, qarr))
    {
        prel = FCLAW2D_PATCH_SAMESIZE;
    }
    else if (q->level < P4EST_QMAXLEVEL &&
             (p4est_quadrant_half_corner_neighbor (q, cornerno, &r),
              p4est_quadrant_exists (p4est, ghost, nt, &r, earr, rparr,
                                     qarr)))
    {
        prel = FCLAW2D_PATCH_HALFSIZE;
    }
    else if (q->level > 0 && p4est_quadrant_child_id (q) == cornerno &&
             (p4est_quadrant_parent (q, &r),
              p4est_quadrant_corner_neighbor (&r, cornerno, &r),
              p4est_quadrant_exists (p4est, ghost, nt, &r, earr, rparr,
                                     qarr)))
    {
        prel = FCLAW2D_PATCH_DOUBLESIZE;
    }
    else
    {
        prel = FCLAW2D_PATCH_BOUNDARY;
    }
    *neighbor_size = prel;

    if (prel != FCLAW2D_PATCH_BOUNDARY)
    {
        P4EST_ASSERT (rparr->elem_count == 1);
        P4EST_ASSERT (qarr->elem_count == 1);
        *rproc = *(int *) sc_array_index (rparr, 0);
        rq = p4est_quadrant_array_index (qarr, 0);
        *rblockno = (p4est_topidx_t) rq->p.piggy3.which_tree;
        *rpatchno = (p4est_topidx_t) rq->p.piggy3.local_num;    /* ghost index */
        P4EST_ASSERT (*rproc == domain->mpirank ||
                      (*rpatchno >= 0
                       && *rpatchno < mesh->ghost_num_quadrants));
        P4EST_ASSERT (*rproc != domain->mpirank
                      || (*rblockno >= 0 && *rblockno < domain->num_blocks
                          && *rpatchno >= 0 && *rpatchno <
                          domain->blocks[*rblockno].num_patches));
    }

    sc_array_reset (earr);
    sc_array_reset (rparr);
    sc_array_reset (qarr);

    return prel != FCLAW2D_PATCH_BOUNDARY;
}

void
fclaw2d_patch_mark_refine (fclaw2d_domain_t * domain, int blockno,
                           int patchno)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    p4est_wrap_mark_refine (wrap,
                            (p4est_locidx_t) blockno,
                            (p4est_locidx_t) patchno);
}

void
fclaw2d_patch_mark_coarsen (fclaw2d_domain_t * domain, int blockno,
                            int patchno)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    p4est_wrap_mark_coarsen (wrap,
                             (p4est_locidx_t) blockno,
                             (p4est_locidx_t) patchno);
}

void
fclaw2d_domain_iterate_adapted (fclaw2d_domain_t * old_domain,
                                fclaw2d_domain_t * new_domain,
                                fclaw2d_match_callback_t mcb, void *user)
{
    int i, oj, nj;
    int oskip, nskip;
    fclaw2d_block_t *old_block, *new_block;
    fclaw2d_patch_t *old_patch, *new_patch;
    fclaw2d_patch_relation_t newsize;

    P4EST_ASSERT (!old_domain->pp_owned);
    P4EST_ASSERT (new_domain->pp_owned);
    P4EST_ASSERT (old_domain->pp == new_domain->pp);
    P4EST_ASSERT (old_domain->num_blocks == new_domain->num_blocks);
    for (i = 0; i < old_domain->num_blocks; i++)
    {
        old_block = old_domain->blocks + i;
        new_block = new_domain->blocks + i;
        for (oj = nj = 0; oj < old_block->num_patches;)
        {
            P4EST_ASSERT (nj < new_block->num_patches);
            old_patch = old_block->patches + oj;
            new_patch = new_block->patches + nj;
            P4EST_ASSERT (abs (old_patch->level - new_patch->level) <= 1);
            if (old_patch->level < new_patch->level)
            {
                /* refinement */
                newsize = FCLAW2D_PATCH_HALFSIZE;
                oskip = 1;
                nskip = P4EST_CHILDREN;
            }
            else if (old_patch->level > new_patch->level)
            {
                /* coarsening */
                newsize = FCLAW2D_PATCH_DOUBLESIZE;
                oskip = P4EST_CHILDREN;
                nskip = 1;
            }
            else
            {
                /* noop */
                newsize = FCLAW2D_PATCH_SAMESIZE;
                oskip = nskip = 1;
            }
            mcb (old_domain, old_patch, new_domain, new_patch,
                 newsize, i, oj, nj, user);

            oj += oskip;
            nj += nskip;
        }
        P4EST_ASSERT (oj == old_block->num_patches);
        P4EST_ASSERT (nj == new_block->num_patches);
    }
}

static void
fclaw2d_domain_assign_for_partition (fclaw2d_domain_t * domain,
                                     void **patch_data)
{
    int blockno, patchno;
    size_t zz;
    fclaw2d_block_t *block;
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_tree_t *tree;
    p4est_quadrant_t *q;

    for (zz = 0, blockno = 0; blockno < domain->num_blocks; ++blockno)
    {
        block = domain->blocks + blockno;
        tree =
            p4est_tree_array_index (wrap->p4est->trees,
                                    (p4est_topidx_t) blockno);

        for (patchno = 0; patchno < block->num_patches; ++zz, ++patchno)
        {
            P4EST_ASSERT (zz ==
                          (size_t) (block->num_patches_before + patchno));

            q = p4est_quadrant_array_index (&tree->quadrants,
                                            (p4est_locidx_t) patchno);
            patch_data[zz] = q->p.user_data;
        }
    }
    P4EST_ASSERT (zz == (size_t) domain->local_num_patches);
}

void
fclaw2d_domain_allocate_before_partition (fclaw2d_domain_t * domain,
                                          size_t data_size,
                                          void ***patch_data)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    P4EST_ASSERT (*patch_data == NULL);

    p4est_reset_data (wrap->p4est, data_size, NULL,
                      wrap->p4est->user_pointer);

    *patch_data = P4EST_ALLOC (void *, domain->local_num_patches);
    fclaw2d_domain_assign_for_partition (domain, *patch_data);
}

void
fclaw2d_domain_retrieve_after_partition (fclaw2d_domain_t * domain,
                                         void ***patch_data)
{
    *patch_data =
        P4EST_REALLOC (*patch_data, void *, domain->local_num_patches);
    fclaw2d_domain_assign_for_partition (domain, *patch_data);
}

void
fclaw2d_domain_free_after_partition (fclaw2d_domain_t * domain,
                                     void ***patch_data)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    P4EST_FREE (*patch_data);
    *patch_data = NULL;

    p4est_reset_data (wrap->p4est, 0, NULL, wrap->p4est->user_pointer);
}

fclaw2d_domain_exchange_t *
fclaw2d_domain_allocate_before_exchange (fclaw2d_domain_t * domain,
                                         size_t data_size)
{
    int i;
    char *m;
    fclaw2d_domain_exchange_t *e;

    e = P4EST_ALLOC (fclaw2d_domain_exchange_t, 1);
    e->data_size = data_size;
    e->num_exchange_patches = domain->num_exchange_patches;
    e->num_ghost_patches = domain->num_ghost_patches;

    e->patch_data = P4EST_ALLOC (void *, domain->num_exchange_patches);
    e->ghost_data = P4EST_ALLOC (void *, domain->num_ghost_patches);
    e->ghost_contiguous_memory = m = P4EST_ALLOC (char,
                                                  (size_t)
                                                  domain->num_ghost_patches *
                                                  data_size);
    for (i = 0; i < domain->num_ghost_patches; ++i)
    {
        e->ghost_data[i] = m;
        m += data_size;
    }

    return e;
}

void
fclaw2d_domain_ghost_exchange (fclaw2d_domain_t * domain,
                               fclaw2d_domain_exchange_t * e)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;

#if 0
    P4EST_LDEBUGF ("Patches exchange %d %d, ghost %d %d\n",
                   e->num_exchange_patches, (int) ghost->mirrors.elem_count,
                   e->num_ghost_patches, (int) ghost->ghosts.elem_count);
#endif

    P4EST_ASSERT (e->num_exchange_patches == (int) ghost->mirrors.elem_count);
    P4EST_ASSERT (e->num_ghost_patches == (int) ghost->ghosts.elem_count);
    p4est_ghost_exchange_custom_data (wrap->p4est, ghost, e->data_size,
                                      e->patch_data,
                                      e->ghost_contiguous_memory);
}

void
fclaw2d_domain_free_after_exchange (fclaw2d_domain_t * domain,
                                    fclaw2d_domain_exchange_t * e)
{
    P4EST_FREE (e->ghost_contiguous_memory);
    P4EST_FREE (e->ghost_data);
    P4EST_FREE (e->patch_data);
    P4EST_FREE (e);
}

void
fclaw2d_domain_serialization_enter (fclaw2d_domain_t * domain)
{
    int mpiret;
    int i;
    MPI_Status status;

    if (domain->mpirank > 0)
    {
        mpiret = MPI_Recv (&i, 1, MPI_INT, domain->mpirank - 1,
                           FCLAW2D_DOMAIN_TAG_SERIALIZE, domain->mpicomm,
                           &status);
        SC_CHECK_MPI (mpiret);
        P4EST_ASSERT (i == 0);
    }
}

void
fclaw2d_domain_serialization_leave (fclaw2d_domain_t * domain)
{
    int mpiret;
    int i = 0;

    if (domain->mpirank + 1 < domain->mpisize)
    {
        mpiret = MPI_Send (&i, 1, MPI_INT, domain->mpirank + 1,
                           FCLAW2D_DOMAIN_TAG_SERIALIZE, domain->mpicomm);
        SC_CHECK_MPI (mpiret);
    }
}
