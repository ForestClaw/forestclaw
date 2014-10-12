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

#include <fclaw2d_convenience.h>
#include <p4est_bits.h>
#include <p4est_wrap.h>

const double fclaw2d_smallest_h = 1. / (double) P4EST_ROOT_LEN;

/** Domain constructor takes ownership of wrap.
 */
static fclaw2d_domain_t *
fclaw2d_domain_new (p4est_wrap_t * wrap, sc_keyvalue_t * attributes)
{
    int mpiret;
    int i, j;
    int level;
    int face;
    int nb, nm, mirror_quadrant_num;
    int local_num_patches;
    int tree_minlevel, local_minlevel;
    int tree_maxlevel, local_maxlevel;
    int levels[2], global_levels[2];
    p4est_qcoord_t qh;
    p4est_connectivity_t *conn = wrap->conn;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
#ifdef FCLAW_ENABLE_DEBUG
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
#endif
    p4est_tree_t *tree;
    p4est_quadrant_t *quad, *mirror;
    fclaw2d_domain_t *domain;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;
    fclaw2d_patch_t *currentbylevel[P4EST_MAXLEVEL + 1];

#ifdef FCLAW_ENABLE_DEBUG
    memset (currentbylevel, 0,
            sizeof (fclaw2d_patch_t *) * (P4EST_MAXLEVEL + 1));
#endif
    domain = FCLAW_ALLOC_ZERO (fclaw2d_domain_t, 1);
    domain->mpicomm = wrap->p4est->mpicomm;
    domain->mpisize = wrap->p4est->mpisize;
    domain->mpirank = wrap->p4est->mpirank;
    domain->pp = (void *) wrap;
    domain->pp_owned = 1;
    domain->attributes = attributes != NULL ? attributes : sc_keyvalue_new ();
    nm = 0;
    domain->num_exchange_patches = (int) ghost->mirrors.elem_count;
    if (domain->num_exchange_patches > 0)
    {
        mirror = p4est_quadrant_array_index (&ghost->mirrors, nm);
        mirror_quadrant_num = (int) mirror->p.piggy3.local_num;
        FCLAW_ASSERT (mirror_quadrant_num >= 0);
        FCLAW_ASSERT (mirror_quadrant_num <
                      (int) wrap->p4est->local_num_quadrants);
    }
    else
    {
        mirror = NULL;
        mirror_quadrant_num = -1;
    }
    domain->num_ghost_patches = (int) ghost->ghosts.elem_count;
    FCLAW_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);
    domain->num_blocks = nb = (int) conn->num_trees;
    domain->blocks = FCLAW_ALLOC_ZERO (fclaw2d_block_t, domain->num_blocks);
    domain->possible_maxlevel = P4EST_QMAXLEVEL;
    local_num_patches = 0;
    local_minlevel = domain->possible_maxlevel;
    local_maxlevel = -1;
    for (i = 0; i < nb; ++i)
    {
        block = domain->blocks + i;
        block->num_patches_before = local_num_patches;
        tree =
            p4est_tree_array_index (wrap->p4est->trees, (p4est_topidx_t) i);
        tree_minlevel = domain->possible_maxlevel;
        tree_maxlevel = -1;
        block->xlower = 0.;
        block->xupper = 1.;
        block->ylower = 0.;
        block->yupper = 1.;
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
            FCLAW_ALLOC_ZERO (fclaw2d_patch_t, block->num_patches);
        block->patchbylevel =
            FCLAW_ALLOC_ZERO (fclaw2d_patch_t *,
                              domain->possible_maxlevel + 1);
        for (j = 0; j < block->num_patches; ++j)
        {
            patch = block->patches + j;
            quad = p4est_quadrant_array_index (&tree->quadrants, (size_t) j);
            patch->level = level = (int) quad->level;
            patch->flags = p4est_quadrant_child_id (quad);
            if (j + P4EST_CHILDREN <= block->num_patches &&
                p4est_quadrant_is_familyv (quad))
            {
                patch->flags |= FCLAW2D_PATCH_FIRST_SIBLING;
            }
            FCLAW_ASSERT (0 <= level && level <= domain->possible_maxlevel);
            if (mirror_quadrant_num == local_num_patches)
            {
                patch->flags |= FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY;
                if (++nm < domain->num_exchange_patches)
                {
                    mirror = p4est_quadrant_array_index (&ghost->mirrors, nm);
                    mirror_quadrant_num = (int) mirror->p.piggy3.local_num;
                    FCLAW_ASSERT (mirror_quadrant_num > local_num_patches);
                    FCLAW_ASSERT (mirror_quadrant_num <
                                  (int) wrap->p4est->local_num_quadrants);
                }
                else
                {
                    mirror = NULL;
                    mirror_quadrant_num = -1;
                }
            }
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
                /* next pointer of previous patch by level in this block */
                FCLAW_ASSERT (currentbylevel[level] != NULL);
                currentbylevel[level]->u.next = patch;
                currentbylevel[level] = patch;
            }
            FCLAW_ASSERT (i == (int) mesh->quad_to_tree[local_num_patches]);
            tree_minlevel = SC_MIN (tree_minlevel, level);
            tree_maxlevel = SC_MAX (tree_maxlevel, level);
            local_num_patches++;
        }
        FCLAW_ASSERT (block->num_patches == 0 ||
                      tree_maxlevel == (int) tree->maxlevel);
        local_minlevel = SC_MIN (local_minlevel, tree_minlevel);
        local_maxlevel = SC_MAX (local_maxlevel, tree_maxlevel);
        block->minlevel = tree_minlevel;
        block->maxlevel = tree_maxlevel;
    }
    FCLAW_ASSERT (local_num_patches ==
                  (int) wrap->p4est->local_num_quadrants);
    FCLAW_ASSERT (nm == domain->num_exchange_patches);
    domain->local_num_patches = local_num_patches;
    domain->local_minlevel = local_minlevel;
    domain->local_maxlevel = local_maxlevel;

    /* allocate ghost patches */
    domain->ghost_patches =
        FCLAW_ALLOC_ZERO (fclaw2d_patch_t, domain->num_ghost_patches);
    for (i = 0; i < domain->num_ghost_patches; ++i)
    {
        patch = domain->ghost_patches + i;
        quad = p4est_quadrant_array_index (&ghost->ghosts, (size_t) i);
        patch->level = level = (int) quad->level;
        patch->flags =
            p4est_quadrant_child_id (quad) | FCLAW2D_PATCH_IS_GHOST;
        FCLAW_ASSERT (0 <= level && level <= domain->possible_maxlevel);
        qh = P4EST_QUADRANT_LEN (level);
        patch->xlower = quad->x * fclaw2d_smallest_h;
        patch->xupper = (quad->x + qh) * fclaw2d_smallest_h;
        patch->ylower = quad->y * fclaw2d_smallest_h;
        patch->yupper = (quad->y + qh) * fclaw2d_smallest_h;
        patch->u.blockno = (int) quad->p.which_tree;
    }

    /* parallel communication of minimum and maximum levels */
    levels[0] = domain->local_minlevel;
    levels[1] = -domain->local_maxlevel;
    mpiret =
        sc_MPI_Allreduce (levels, global_levels, 2, sc_MPI_INT, sc_MPI_MIN,
                          domain->mpicomm);
    SC_CHECK_MPI (mpiret);
    domain->global_minlevel = global_levels[0];
    domain->global_maxlevel = -global_levels[1];
    FCLAW_ASSERT (0 <= domain->global_minlevel);
    FCLAW_ASSERT (domain->global_minlevel <= domain->global_maxlevel);
    FCLAW_ASSERT (domain->global_maxlevel <= domain->possible_maxlevel);
    domain->global_num_patches = (int64_t) wrap->p4est->global_num_quadrants;
    domain->global_num_patches_before =
        (int64_t) wrap->p4est->global_first_quadrant[domain->mpirank];

    return domain;
}

static void
fclaw2d_check_initial_level (sc_MPI_Comm mpicomm, int initial_level)
{
    int mpiret;
    int rank;

    mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORTF (rank != 0 || initial_level <= P4EST_QMAXLEVEL,
                     "Initial level %d too fine for p4est", initial_level);
}

fclaw2d_domain_t *
fclaw2d_domain_new_unitsquare (sc_MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return fclaw2d_domain_new (p4est_wrap_new_unitsquare (mpicomm,
                                                          initial_level),
                               NULL);
}

fclaw2d_domain_t *
fclaw2d_domain_new_torus (sc_MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return
        fclaw2d_domain_new (p4est_wrap_new_periodic (mpicomm, initial_level),
                            NULL);
}

fclaw2d_domain_t *
fclaw2d_domain_new_twosphere (sc_MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return
        fclaw2d_domain_new (p4est_wrap_new_pillow (mpicomm, initial_level),
                            NULL);
}

fclaw2d_domain_t *
fclaw2d_domain_new_cubedsphere (sc_MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return fclaw2d_domain_new (p4est_wrap_new_cubed (mpicomm, initial_level),
                               NULL);
}

fclaw2d_domain_t *
fclaw2d_domain_new_disk (sc_MPI_Comm mpicomm, int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return fclaw2d_domain_new (p4est_wrap_new_disk (mpicomm, initial_level),
                               NULL);
}

fclaw2d_domain_t *
fclaw2d_domain_new_brick_map (sc_MPI_Comm mpicomm,
                              int blocks_in_x, int blocks_in_y,
                              int periodic_in_x, int periodic_in_y,
                              int initial_level, fclaw2d_map_context_t * cont)
{
    p4est_wrap_t *wrap;
    fclaw2d_domain_t *domain;

    fclaw2d_check_initial_level (mpicomm, initial_level);
    wrap = p4est_wrap_new_brick (mpicomm, blocks_in_x, blocks_in_y,
                                 periodic_in_x, periodic_in_y, initial_level);
    domain = fclaw2d_domain_new (wrap, NULL);
    if (cont != NULL)
    {
        fclaw2d_domain_attribute_add (domain, "fclaw_map_context", cont);
    }

    return domain;
}

fclaw2d_domain_t *
fclaw2d_domain_new_conn_map (sc_MPI_Comm mpicomm, int initial_level,
                             p4est_connectivity_t * conn,
                             fclaw2d_map_context_t * cont)
{
    p4est_wrap_t *wrap;
    fclaw2d_domain_t *domain;

    fclaw2d_check_initial_level (mpicomm, initial_level);
    wrap = p4est_wrap_new_conn (mpicomm, conn, initial_level);
    domain = fclaw2d_domain_new (wrap, NULL);
    fclaw2d_domain_attribute_add (domain, "fclaw_map_context", cont);

    return domain;
}

void
fclaw2d_domain_destroy (fclaw2d_domain_t * domain)
{
    int i;
    fclaw2d_block_t *block;

    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        FCLAW_FREE (block->patches);
        FCLAW_FREE (block->patchbylevel);
    }
    FCLAW_FREE (domain->blocks);

    FCLAW_FREE (domain->ghost_patches);

    if (domain->pp_owned)
    {
        p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
        p4est_wrap_destroy (wrap);
        sc_keyvalue_destroy (domain->attributes);
    }
    FCLAW_FREE (domain);
}

fclaw2d_domain_t *
fclaw2d_domain_adapt (fclaw2d_domain_t * domain)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (domain->pp_owned);
    if (p4est_wrap_adapt (wrap))
    {
        domain->pp_owned = 0;
        return fclaw2d_domain_new (wrap, domain->attributes);
    }
    else
    {
        return NULL;
    }
}

fclaw2d_domain_t *
fclaw2d_domain_partition (fclaw2d_domain_t * domain, int weight_exponent)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (domain->pp_owned);
    if (p4est_wrap_partition (wrap, weight_exponent))
    {
        domain->pp_owned = 0;
        return fclaw2d_domain_new (wrap, domain->attributes);
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

    FCLAW_ASSERT (domain->pp_owned);
    p4est_wrap_complete (wrap);
}

static void
fclaw2d_domain_list_level_callback (fclaw2d_domain_t * domain,
                                    fclaw2d_patch_t * patch, int block_no,
                                    int patch_no, void *user)
{
    FCLAW_ASSERT (0 <= block_no && block_no < domain->num_blocks);
    FCLAW_ASSERT (0 <= patch_no &&
                  patch_no < domain->blocks[block_no].num_patches);
    FCLAW_ASSERT (patch == domain->blocks[block_no].patches + patch_no);

    (*(int *) user)++;
}

void
fclaw2d_domain_list_levels (fclaw2d_domain_t * domain, int lp)
{
    int level;
    int count, count_all;

    P4EST_LOGF (lp, "Local minimum/maximum levels: %2d %2d\n",
                domain->local_minlevel, domain->local_maxlevel);
    P4EST_GLOBAL_LOGF (lp, "Global minimum/maximum levels: %2d %2d\n",
                       domain->global_minlevel, domain->global_maxlevel);
    count_all = 0;
    for (level = domain->local_minlevel; level <= domain->local_maxlevel;
         ++level)
    {
        count = 0;
        fclaw2d_domain_iterate_level (domain, level,
                                      fclaw2d_domain_list_level_callback,
                                      &count);
        P4EST_LOGF (lp, "Patches on level %2d: %9d\n", level, count);
        count_all += count;
    }
    FCLAW_ASSERT (count_all == domain->local_num_patches);
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
    int faceno, cornerno, rcorner;
    int rproc[2], rblockno, rpatchno[2], rfaceno;

    FCLAW_ASSERT (0 <= block_no && block_no < domain->num_blocks);
    FCLAW_ASSERT (0 <= patch_no &&
                  patch_no < domain->blocks[block_no].num_patches);
    FCLAW_ASSERT (patch == domain->blocks[block_no].patches + patch_no);

    for (faceno = 0; faceno < P4EST_FACES; ++faceno)
    {
        fnt = fclaw2d_patch_face_neighbors (domain, block_no, patch_no,
                                            faceno, rproc, &rblockno,
                                            rpatchno, &rfaceno);
        P4EST_LOGF (ln->lp, "Block %d patch %d face %d neighbor %d\n",
                    block_no, patch_no, faceno, (int) fnt);
    }
    for (cornerno = 0; cornerno < P4EST_CHILDREN; ++cornerno)
    {
        (void) fclaw2d_patch_corner_neighbors (domain, block_no, patch_no,
                                               cornerno, rproc, &rblockno,
                                               rpatchno, &rcorner, &fnt);
        P4EST_LOGF (ln->lp, "Block %d patch %d corner %d %d neighbor %d\n",
                    block_no, patch_no, cornerno, rcorner, (int) fnt);
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
    FCLAW_ASSERT (ln.count == domain->local_num_patches);
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
    int k;

    FCLAW_ASSERT (old_domain->pp == new_domain->pp);
    FCLAW_ASSERT (old_domain->num_blocks == new_domain->num_blocks);
    FCLAW_ASSERT (0 <= blockno && blockno < old_domain->num_blocks);
    FCLAW_ASSERT (0 <= old_patchno &&
                  old_patchno < old_domain->blocks[blockno].num_patches);
    FCLAW_ASSERT (0 <= new_patchno &&
                  new_patchno < new_domain->blocks[blockno].num_patches);
    FCLAW_ASSERT (old_patch ==
                  old_domain->blocks[blockno].patches + old_patchno);
    FCLAW_ASSERT (new_patch ==
                  new_domain->blocks[blockno].patches + new_patchno);

    if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        /* refinement */
        FCLAW_ASSERT (new_patchno + P4EST_CHILDREN <=
                      new_domain->blocks[blockno].num_patches);
        for (k = 0; k < P4EST_CHILDREN; ++k)
        {
            FCLAW_ASSERT (new_patch[k].level == old_patch->level + 1);
            FCLAW_ASSERT (fclaw2d_patch_childid (&new_patch[k]) == k);
        }
        P4EST_LOGF (lp, "Block %d refinement %d to %d\n",
                    blockno, old_patchno, new_patchno);
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        /* coarsening */
        FCLAW_ASSERT (old_patchno + P4EST_CHILDREN <=
                      old_domain->blocks[blockno].num_patches);
        for (k = 0; k < P4EST_CHILDREN; ++k)
        {
            FCLAW_ASSERT (old_patch[k].level == new_patch->level + 1);
            FCLAW_ASSERT (fclaw2d_patch_childid (&old_patch[k]) == k);
        }
        P4EST_LOGF (lp, "Block %d coarsening %d to %d\n",
                    blockno, old_patchno, new_patchno);
    }
    else
    {
        /* noop */
        FCLAW_ASSERT (newsize == FCLAW2D_PATCH_SAMESIZE);
        FCLAW_ASSERT (old_patch->level == new_patch->level);
        FCLAW_ASSERT (fclaw2d_patch_childid (old_patch) ==
                      fclaw2d_patch_childid (new_patch));
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
