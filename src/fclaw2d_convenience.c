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

#include "fclaw2d_convenience.h"
#include <p4est_bits.h>
#include <p4est_wrap.h>

const double fclaw2d_root_len = (double) P4EST_ROOT_LEN;
const double fclaw2d_smallest_h = 1. / (double) P4EST_ROOT_LEN;

static void
fclaw2d_domain_mcp (const double xyc[2], double xyzp[P4EST_DIM],
                    fclaw2d_domain_t * domain, void *user)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_connectivity_t *conn = wrap->conn;
    p4est_topidx_t treeid;

    treeid = (p4est_topidx_t) (long) user;
    P4EST_ASSERT (0 <= treeid && treeid < conn->num_trees);

    p4est_qcoord_to_vertex (conn, treeid,
                            (p4est_qcoord_t) (xyc[0] * fclaw2d_root_len),
                            (p4est_qcoord_t) (xyc[1] * fclaw2d_root_len),
                            xyzp);
}

/** Domain constructor takes ownership of wrap.
 */
static fclaw2d_domain_t *
fclaw2d_domain_new (p4est_wrap_t * wrap)
{
    int mpiret;
    int i, j;
    int level;
    int face;
    int nb;
    int num_patches_all;
    int tree_minlevel, minlevel_all;
    int tree_maxlevel, maxlevel_all;
    int levels[2], global_levels[2];
    p4est_qcoord_t qh;
    p4est_connectivity_t *conn = wrap->conn;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    fclaw2d_domain_t *domain;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;
    fclaw2d_patch_t *currentbylevel[P4EST_MAXLEVEL + 1];

#ifdef P4EST_DEBUG
    memset (currentbylevel, 0,
            sizeof (fclaw2d_patch_t *) * (P4EST_MAXLEVEL + 1));
#endif
    domain = P4EST_ALLOC_ZERO (fclaw2d_domain_t, 1);
    domain->mpicomm = wrap->p4est->mpicomm;
    domain->mpisize = wrap->p4est->mpisize;
    domain->mpirank = wrap->p4est->mpirank;
    domain->pp = (void *) wrap;
    domain->pp_owned = 1;
    domain->num_blocks = nb = (int) conn->num_trees;
    domain->blocks = P4EST_ALLOC_ZERO (fclaw2d_block_t, domain->num_blocks);
    domain->patch_to_block =
        P4EST_ALLOC (int, wrap->p4est->local_num_quadrants);
    domain->possible_maxlevel = P4EST_QMAXLEVEL;
    num_patches_all = 0;
    minlevel_all = domain->possible_maxlevel;
    maxlevel_all = 0;
    for (i = 0; i < nb; ++i)
    {
        block = domain->blocks + i;
        block->num_patches_before = num_patches_all;
        tree =
            p4est_tree_array_index (wrap->p4est->trees, (p4est_topidx_t) i);
        tree_minlevel = domain->possible_maxlevel;
        tree_maxlevel = 0;
        block->xlower = 0.;
        block->xupper = 1.;
        block->ylower = 0.;
        block->yupper = 1.;
        if (conn->num_vertices > 0)
        {
            P4EST_ASSERT (conn->vertices != NULL);
            block->mapc2m = fclaw2d_domain_mcp;
            block->mapc2m_user = (void *) (long) i;
        }
        for (face = 0; face < P4EST_FACES; ++face)
        {
            if (conn->tree_to_tree[P4EST_FACES * i + face] ==
                (p4est_topidx_t) i
                && conn->tree_to_face[P4EST_FACES * i + face] ==
                (int8_t) face)
            {
                block->is_boundary[face] = 1;
            }
        }
        block->num_patches = (int) tree->quadrants.elem_count;
        block->patches =
            P4EST_ALLOC_ZERO (fclaw2d_patch_t, block->num_patches);
        block->patchbylevel =
            P4EST_ALLOC_ZERO (fclaw2d_patch_t *,
                              domain->possible_maxlevel + 1);
        for (j = 0; j < block->num_patches; ++j)
        {
            patch = block->patches + j;
            quad = p4est_quadrant_array_index (&tree->quadrants, (size_t) j);
            patch->level = level = (int) quad->level;
            patch->childid = p4est_quadrant_child_id (quad);
            P4EST_ASSERT (0 <= level && level <= domain->possible_maxlevel);
            qh = P4EST_QUADRANT_LEN (level);
            patch->xlower = quad->x * fclaw2d_smallest_h;
            patch->xupper = (quad->x + qh) * fclaw2d_smallest_h;
            patch->ylower = quad->y * fclaw2d_smallest_h;
            patch->yupper = (quad->y + qh) * fclaw2d_smallest_h;
            if (block->patchbylevel[level] == NULL)
            {
                /* this is the first patch of this level in this block */
                block->patchbylevel[level] = currentbylevel[level] = patch;
            }
            else
            {
                /* assign next pointer of previous patch by level in this block */
                P4EST_ASSERT (currentbylevel[level] != NULL);
                currentbylevel[level]->next = patch;
                currentbylevel[level] = patch;
            }
            domain->patch_to_block[num_patches_all++] = i;
            tree_minlevel = SC_MIN (tree_minlevel, level);
            tree_maxlevel = SC_MAX (tree_maxlevel, level);
        }
        P4EST_ASSERT (tree_maxlevel == (int) tree->maxlevel);
        minlevel_all = SC_MIN (minlevel_all, tree_minlevel);
        maxlevel_all = SC_MAX (maxlevel_all, tree_maxlevel);
        block->minlevel = tree_minlevel;
        block->maxlevel = tree_maxlevel;
    }
    P4EST_ASSERT (num_patches_all == (int) wrap->p4est->local_num_quadrants);
    domain->num_patches_all = num_patches_all;
    domain->minlevel_all = minlevel_all;
    domain->maxlevel_all = maxlevel_all;

    /* parallel communication of minimum and maximum levels */
    levels[0] = domain->minlevel_all;
    levels[1] = -domain->maxlevel_all;
    mpiret = MPI_Allreduce (levels, global_levels, 2, MPI_INT, MPI_MIN,
                            domain->mpicomm);
    SC_CHECK_MPI (mpiret);
    domain->global_minlevel = global_levels[0];
    domain->global_maxlevel = -global_levels[1];

    return domain;
}

static void
fclaw2d_check_initial_level (MPI_Comm mpicomm, int initial_level)
{
    int mpiret;
    int rank;

    mpiret = MPI_Comm_rank (mpicomm, &rank);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORTF (rank != 0 || initial_level <= P4EST_QMAXLEVEL,
                     "Initial level %d too fine for p4est", initial_level);
}

fclaw2d_domain_t *
fclaw2d_domain_new_unitsquare (MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return fclaw2d_domain_new (p4est_wrap_new_unitsquare (mpicomm,
                                                          initial_level));
}

fclaw2d_domain_t *
fclaw2d_domain_new_pillow (MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return
        fclaw2d_domain_new (p4est_wrap_new_pillow (mpicomm, initial_level));
}

fclaw2d_domain_t *
fclaw2d_domain_new_moebius (MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return
        fclaw2d_domain_new (p4est_wrap_new_moebius (mpicomm, initial_level));
}

void
fclaw2d_domain_destroy (fclaw2d_domain_t * domain)
{
    int i;
    fclaw2d_block_t *block;

    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        P4EST_FREE (block->patches);
        P4EST_FREE (block->patchbylevel);
    }
    P4EST_FREE (domain->patch_to_block);
    P4EST_FREE (domain->blocks);

    if (domain->pp_owned)
    {
        p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
        p4est_wrap_destroy (wrap);
    }
    P4EST_FREE (domain);
}

fclaw2d_domain_t *
fclaw2d_domain_adapt (fclaw2d_domain_t * domain)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    P4EST_ASSERT (domain->pp_owned);
    if (p4est_wrap_refine (wrap))
    {
        domain->pp_owned = 0;
        return fclaw2d_domain_new (wrap);
    }
    else
    {
        return NULL;
    }
}

fclaw2d_domain_t *
fclaw2d_domain_partition (fclaw2d_domain_t * domain)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    P4EST_ASSERT (domain->pp_owned);
    if (p4est_wrap_partition (wrap))
    {
        domain->pp_owned = 0;
        return fclaw2d_domain_new (wrap);
    }
    else
    {
        return NULL;
    }
}

void
fclaw2d_domain_complete (fclaw2d_domain_t * domain)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    P4EST_ASSERT (domain->pp_owned);
    p4est_wrap_complete (wrap);
}

static void
fclaw2d_domain_list_level_callback (fclaw2d_domain_t * domain,
                                    fclaw2d_patch_t * patch, int block_no,
                                    int patch_no, void *user)
{
    P4EST_ASSERT (0 <= block_no && block_no < domain->num_blocks);
    P4EST_ASSERT (0 <= patch_no &&
                  patch_no < domain->blocks[block_no].num_patches);
    P4EST_ASSERT (patch == domain->blocks[block_no].patches + patch_no);

    (*(int *) user)++;
}

void
fclaw2d_domain_list_levels (fclaw2d_domain_t * domain, int lp)
{
    int level;
    int count, count_all;

    P4EST_ASSERT (0 <= domain->minlevel_all &&
                  domain->minlevel_all <= domain->possible_maxlevel);
    P4EST_ASSERT (0 <= domain->maxlevel_all &&
                  domain->maxlevel_all <= domain->possible_maxlevel);
    P4EST_ASSERT (domain->minlevel_all <= domain->maxlevel_all);
    P4EST_GLOBAL_LOGF (lp, "Global minimum/maximum levels: %2d %2d\n",
                       domain->global_minlevel, domain->global_maxlevel);
    P4EST_LOGF (lp, "Minimum/maximum levels: %2d %2d\n",
                domain->minlevel_all, domain->maxlevel_all);
    count_all = 0;
    for (level = domain->minlevel_all; level <= domain->maxlevel_all; ++level)
    {
        count = 0;
        fclaw2d_domain_iterate_level (domain, level,
                                      fclaw2d_domain_list_level_callback,
                                      &count);
        P4EST_LOGF (lp, "Patches on level %2d: %9d\n", level, count);
        count_all += count;
    }
    P4EST_ASSERT (count_all == domain->num_patches_all);
}

typedef struct fclaw2d_domain_list_neighbors
{
    int lp;
    int count;
}
fclaw2d_domain_list_neighbors_t;

static void
fclaw2d_domain_list_neighbors_callback (fclaw2d_domain_t * domain,
                                        fclaw2d_patch_t * patch, int block_no,
                                        int patch_no, void *user)
{
    fclaw2d_domain_list_neighbors_t *ln = user;
    fclaw2d_patch_relation_t fnt;
    int hasnb;
    int faceno, cornerno;
    int rproc[2], rblockno, rpatchno[2], rfaceno;

    P4EST_ASSERT (0 <= block_no && block_no < domain->num_blocks);
    P4EST_ASSERT (0 <= patch_no &&
                  patch_no < domain->blocks[block_no].num_patches);
    P4EST_ASSERT (patch == domain->blocks[block_no].patches + patch_no);

    for (faceno = 0; faceno < P4EST_FACES; ++faceno)
    {
        fnt = fclaw2d_patch_face_neighbors (domain, block_no, patch_no,
                                            faceno, rproc, &rblockno,
                                            rpatchno, &rfaceno);
        P4EST_LOGF (ln->lp, "Block %d patch %d face %d neighbor %d\n",
                    block_no, patch_no, faceno, fnt);
    }
    for (cornerno = 0; cornerno < P4EST_CHILDREN; ++cornerno)
    {
        hasnb = fclaw2d_patch_corner_neighbors (domain, block_no, patch_no,
                                                cornerno, rproc, &rblockno,
                                                rpatchno, &fnt);
        P4EST_LOGF (ln->lp, "Block %d patch %d corner %d neighbor %d\n",
                    block_no, patch_no, cornerno, fnt);
    }

    ++ln->count;
}

void
fclaw2d_domain_list_neighbors (fclaw2d_domain_t * domain, int lp)
{
    fclaw2d_domain_list_neighbors_t ln;

    ln.lp = lp;
    ln.count = 0;
    fclaw2d_domain_iterate_patches (domain,
                                    fclaw2d_domain_list_neighbors_callback,
                                    &ln);
    P4EST_ASSERT (ln.count == domain->num_patches_all);
}

static void
fclaw2d_domain_list_adapted_callback (fclaw2d_domain_t * old_domain,
                                      fclaw2d_patch_t * old_patch,
                                      fclaw2d_domain_t * new_domain,
                                      fclaw2d_patch_t * new_patch,
                                      fclaw2d_patch_relation_t newsize,
                                      int blockno,
                                      int old_patchno, int new_patchno,
                                      void *user)
{
    int lp = *(int *) user;
    int level;
    int k;

    P4EST_ASSERT (old_domain->pp == new_domain->pp);
    P4EST_ASSERT (old_domain->num_blocks == new_domain->num_blocks);
    P4EST_ASSERT (0 <= blockno && blockno < old_domain->num_blocks);
    P4EST_ASSERT (0 <= old_patchno &&
                  old_patchno < old_domain->blocks[blockno].num_patches);
    P4EST_ASSERT (0 <= new_patchno &&
                  new_patchno < new_domain->blocks[blockno].num_patches);
    P4EST_ASSERT (old_patch ==
                  old_domain->blocks[blockno].patches + old_patchno);
    P4EST_ASSERT (new_patch ==
                  new_domain->blocks[blockno].patches + new_patchno);

    if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        /* refinement */
        level = old_patch->level + 1;
        P4EST_ASSERT (new_patchno + P4EST_CHILDREN <=
                      new_domain->blocks[blockno].num_patches);
        for (k = 0; k < P4EST_CHILDREN; ++k)
        {
            P4EST_ASSERT (new_patch[k].level == level);
            P4EST_ASSERT (new_patch[k].childid == k);
        }
        P4EST_LOGF (lp, "Block %d refinement %d to %d\n",
                    blockno, old_patchno, new_patchno);
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        /* coarsening */
        level = new_patch->level + 1;
        P4EST_ASSERT (old_patchno + P4EST_CHILDREN <=
                      old_domain->blocks[blockno].num_patches);
        for (k = 0; k < P4EST_CHILDREN; ++k)
        {
            P4EST_ASSERT (old_patch[k].level == level);
            P4EST_ASSERT (old_patch[k].childid == k);
        }
        P4EST_LOGF (lp, "Block %d coarsening %d to %d\n",
                    blockno, old_patchno, new_patchno);
    }
    else
    {
        /* noop */
        P4EST_ASSERT (newsize == FCLAW2D_PATCH_SAMESIZE);
        P4EST_ASSERT (old_patch->level == new_patch->level);
        P4EST_ASSERT (old_patch->childid == new_patch->childid);
        P4EST_LOGF (lp, "Block %d noop patch %d and %d\n",
                    blockno, old_patchno, new_patchno);
    }
}

void
fclaw2d_domain_list_adapted (fclaw2d_domain_t * old_domain,
                             fclaw2d_domain_t * new_domain, int log_priority)
{
    fclaw2d_domain_iterate_adapted (old_domain, new_domain,
                                    fclaw2d_domain_list_adapted_callback,
                                    (void *) &log_priority);
}
