/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#include <sc_notify.h>
#ifndef P4_TO_P8
#include <fclaw2d_convenience.h>
#include <p4est_bits.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_wrap.h>
#else
#include <fclaw3d_convenience.h>
#include <p8est_bits.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_wrap.h>

// for dimension dependent values
// needs to be defined AFTER all other headers
#define d2 d3
#define fclaw_patch_bounds_2d_t fclaw_patch_bounds_3d_t
#define fclaw_block_d2_t fclaw_block_d3_t

#endif

const double fclaw2d_smallest_h = 1. / (double) P4EST_ROOT_LEN;

static void
fclaw2d_patch_set_boundary_xylower (fclaw_patch_t * patch,
                                    p4est_quadrant_t * quad)
{
#ifndef P4_TO_P8
    patch->dim = 2;
#else
    patch->dim = 3;
#endif

    p4est_qcoord_t qh;

    qh = P4EST_QUADRANT_LEN (quad->level);
    if (quad->x == 0)
    {
        patch->flags |= FCLAW2D_PATCH_ON_BLOCK_FACE_0;
    }
    if (quad->x + qh == P4EST_ROOT_LEN)
    {
        patch->flags |= FCLAW2D_PATCH_ON_BLOCK_FACE_1;
    }
    if (quad->y == 0)
    {
        patch->flags |= FCLAW2D_PATCH_ON_BLOCK_FACE_2;
    }
    if (quad->y + qh == P4EST_ROOT_LEN)
    {
        patch->flags |= FCLAW2D_PATCH_ON_BLOCK_FACE_3;
    }
#ifdef P4_TO_P8
    if (quad->z == 0)
    {
        patch->flags |= FCLAW3D_PATCH_ON_BLOCK_FACE_4;
    }
    if (quad->z + qh == P4EST_ROOT_LEN)
    {
        patch->flags |= FCLAW3D_PATCH_ON_BLOCK_FACE_5;
    }
#endif
    /* This suffices to test for block corners by using bitwise and */

    patch->d2->xlower = quad->x * fclaw2d_smallest_h;
    patch->d2->xupper = (quad->x + qh) * fclaw2d_smallest_h;
    patch->d2->ylower = quad->y * fclaw2d_smallest_h;
    patch->d2->yupper = (quad->y + qh) * fclaw2d_smallest_h;
#ifdef P4_TO_P8
    patch->d3->zlower = quad->z * fclaw2d_smallest_h;
    patch->d3->zupper = (quad->z + qh) * fclaw2d_smallest_h;
#endif
}

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
    int block_nm_pre;
    int local_num_patches;
    int tree_minlevel, local_minlevel;
    int tree_maxlevel, local_maxlevel;
    int levels[2], global_levels[2];
    p4est_topidx_t vnum;
    p4est_connectivity_t *conn = wrap->conn;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
#ifdef FCLAW_ENABLE_DEBUG
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
#endif
    p4est_tree_t *tree;
    p4est_quadrant_t *quad, *mirror;
    fclaw2d_domain_t *domain;
    fclaw2d_block_t *block;
    fclaw_patch_t *patch;
    fclaw_patch_t *currentbylevel[P4EST_MAXLEVEL + 1];

#ifdef FCLAW_ENABLE_DEBUG
    memset (currentbylevel, 0,
            sizeof (fclaw_patch_t *) * (P4EST_MAXLEVEL + 1));
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
        domain->exchange_patches = FCLAW_ALLOC (fclaw_patch_t *,
                                                domain->num_exchange_patches);
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

    /* prepare propagation of refinement/coarsening marks */
    domain->p.smooth_refine = 0;
    domain->p.smooth_level = 0;
    domain->mirror_target_levels =
        FCLAW_ALLOC (void *, domain->num_exchange_patches);
    domain->ghost_target_levels =
        FCLAW_ALLOC (int, domain->num_ghost_patches);

    /* loop through all local patches, block by block */
    for (i = 0; i < nb; ++i)
    {
        block = domain->blocks + i;
        block->num_patches_before = local_num_patches;
        tree =
            p4est_tree_array_index (wrap->p4est->trees, (p4est_topidx_t) i);
        tree_minlevel = domain->possible_maxlevel;
        tree_maxlevel = -1;
        block->d2 = FCLAW_ALLOC_ZERO (fclaw_block_d2_t, 1);
        block->d2->xlower = 0.;
        block->d2->xupper = 1.;
        block->d2->ylower = 0.;
        block->d2->yupper = 1.;
#ifdef P4_TO_P8
        block->d3->zlower = 0.;
        block->d3->zupper = 1.;
#endif
        if (conn->vertices != NULL && conn->tree_to_vertex != NULL)
        {
            for (j = 0; j < P4EST_CHILDREN; ++j)
            {
                vnum = conn->tree_to_vertex[P4EST_CHILDREN * i + j];
                FCLAW_ASSERT (0 <= vnum && vnum < conn->num_vertices);
                memcpy (block->d2->vertices + 3 * j, conn->vertices + 3 * vnum,
                        3 * sizeof (double));
            }
        }
        else
        {
            memset (block->d2->vertices, 0, P4EST_CHILDREN * 3 * sizeof (double));
        }
        for (face = 0; face < P4EST_FACES; ++face)
        {
            if (conn->tree_to_tree[P4EST_FACES * i + face] ==
                (p4est_topidx_t) i
                && conn->tree_to_face[P4EST_FACES * i + face] ==
                (int8_t) face)
            {
                block->d2->is_boundary[face] = 1;
            }
        }
        block->num_patches = (int) tree->quadrants.elem_count;
        block->patches =
            FCLAW_ALLOC_ZERO (fclaw_patch_t, block->num_patches);
        block->patch_bounds =
            FCLAW_ALLOC_ZERO (fclaw_patch_bounds_2d_t, block->num_patches);
        block->patchbylevel =
            FCLAW_ALLOC_ZERO (fclaw_patch_t *,
                              domain->possible_maxlevel + 1);

        block_nm_pre = nm;
        for (j = 0; j < block->num_patches; ++j)
        {
            patch = block->patches + j;
            patch->d2 = block->patch_bounds + j;
            quad = p4est_quadrant_array_index (&tree->quadrants, (size_t) j);
            patch->target_level = patch->level = level = (int) quad->level;
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
                domain->exchange_patches[nm] = patch;
                domain->mirror_target_levels[nm] = &patch->target_level;
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
            fclaw2d_patch_set_boundary_xylower (patch, quad);
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
        block->num_exchange_patches = nm - block_nm_pre;
        if (block->num_exchange_patches > 0)
        {
            block->exchange_patches = domain->exchange_patches + block_nm_pre;
        }
    }
    FCLAW_ASSERT (local_num_patches ==
                  (int) wrap->p4est->local_num_quadrants);
    FCLAW_ASSERT (nm == domain->num_exchange_patches);
    domain->local_num_patches = local_num_patches;
    domain->local_minlevel = local_minlevel;
    domain->local_maxlevel = local_maxlevel;

    /* allocate ghost patches */
    domain->ghost_patches =
        FCLAW_ALLOC_ZERO (fclaw_patch_t, domain->num_ghost_patches);
    domain->ghost_patch_bounds =
        FCLAW_ALLOC_ZERO (fclaw_patch_bounds_2d_t,
                          domain->num_ghost_patches);
    for (i = 0; i < domain->num_ghost_patches; ++i)
    {
        patch = domain->ghost_patches + i;
        patch->d2 = domain->ghost_patch_bounds + i;
        quad = p4est_quadrant_array_index (&ghost->ghosts, (size_t) i);
        patch->target_level = patch->level = level = (int) quad->level;
        patch->flags =
            p4est_quadrant_child_id (quad) | FCLAW2D_PATCH_IS_GHOST;
        FCLAW_ASSERT (0 <= level && level <= domain->possible_maxlevel);
        fclaw2d_patch_set_boundary_xylower (patch, quad);
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

#ifndef P4_TO_P8

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
fclaw2d_domain_new_disk (sc_MPI_Comm mpicomm,
                         int periodic_in_x, int periodic_in_y,
                         int initial_level)
{
    fclaw2d_check_initial_level (mpicomm, initial_level);
    return fclaw2d_domain_new
        (p4est_wrap_new_disk (mpicomm, periodic_in_x, periodic_in_y,
                              initial_level), NULL);
}

#endif /* P4_TO_P8 */

fclaw2d_domain_t *
fclaw2d_domain_new_brick (sc_MPI_Comm mpicomm,
                          int blocks_in_x, int blocks_in_y,
#ifdef P4_TO_P8
                          int blocks_in_z,
#endif
                          int periodic_in_x, int periodic_in_y,
#ifdef P4_TO_P8
                          int periodic_in_z,
#endif
                          int initial_level)
{
    p4est_wrap_t *wrap;

    fclaw2d_check_initial_level (mpicomm, initial_level);
    wrap = p4est_wrap_new_brick (mpicomm, blocks_in_x, blocks_in_y,
#ifdef P4_TO_P8
                                 blocks_in_z,
#endif
                                 periodic_in_x, periodic_in_y,
#ifdef P4_TO_P8
                                 periodic_in_z,
#endif
                                 initial_level);
    return fclaw2d_domain_new (wrap, NULL);
}

fclaw2d_domain_t *
fclaw2d_domain_new_conn (sc_MPI_Comm mpicomm, int initial_level,
                         p4est_connectivity_t * conn)
{
    p4est_wrap_t *wrap;
    fclaw2d_domain_t *domain;

    fclaw2d_check_initial_level (mpicomm, initial_level);
    wrap = p4est_wrap_new_conn (mpicomm, conn, initial_level);
    domain = fclaw2d_domain_new (wrap, NULL);

    return domain;
}

#ifndef P4_TO_P8

/* function to be removed once no longer called by applications */

fclaw2d_domain_t *
fclaw2d_domain_new_conn_map (sc_MPI_Comm mpicomm, int initial_level,
                             p4est_connectivity_t * conn,
                             fclaw2d_map_context_t * cont)
{
    fclaw2d_domain_t *domain =
      fclaw2d_domain_new_conn (mpicomm, initial_level, conn);

    fclaw2d_domain_attribute_add (domain, "fclaw_map_context", cont);

    return domain;
}

#endif

void
fclaw2d_domain_destroy (fclaw2d_domain_t * domain)
{
    int i;
    fclaw2d_block_t *block;

    FCLAW_ASSERT (!domain->just_adapted);
    FCLAW_ASSERT (!domain->just_partitioned);

    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        FCLAW_FREE (block->patches);
        FCLAW_FREE (block->patch_bounds);
        FCLAW_FREE (block->patchbylevel);
    }
    FCLAW_FREE (domain->blocks);

    FCLAW_FREE (domain->ghost_patches);
    FCLAW_FREE (domain->ghost_patch_bounds);
    FCLAW_FREE (domain->ghost_target_levels);
    FCLAW_FREE (domain->mirror_target_levels);
    FCLAW_FREE (domain->exchange_patches);

    if (domain->pp_owned)
    {
        p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
        p4est_wrap_destroy (wrap);
        sc_keyvalue_destroy (domain->attributes);
    }
    FCLAW_FREE (domain);
}

static void
fclaw2d_domain_copy_parameters (fclaw2d_domain_t * domain_dest,
                                fclaw2d_domain_t * domain_src)
{
    memcpy (&domain_dest->p, &domain_src->p,
            sizeof (fclaw2d_domain_persist_t));
}

static fclaw_patch_t *
fclaw2d_domain_get_neighbor_patch (fclaw2d_domain_t * domain,
                                   int nproc, int nblockno, int npatchno)
{
    fclaw2d_block_t *block;

    /* the block number should always be right */
    FCLAW_ASSERT (0 <= nblockno && nblockno < domain->num_blocks);

    /* is this a ghost patch */
    if (nproc != domain->mpirank)
    {
        FCLAW_ASSERT (0 <= npatchno && npatchno < domain->num_ghost_patches);
        FCLAW_ASSERT (domain->ghost_patches[npatchno].u.blockno == nblockno);
        return domain->ghost_patches + npatchno;
    }

    /* local patch */
    block = domain->blocks + nblockno;
    FCLAW_ASSERT (0 <= npatchno && npatchno < block->num_patches);
    return block->patches + npatchno;
}

fclaw2d_domain_t *
fclaw2d_domain_adapt (fclaw2d_domain_t * domain)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (domain->pp_owned);
    P4EST_ASSERT (!domain->just_adapted);
    P4EST_ASSERT (!domain->just_partitioned);

    /* propagate desired refinement level to neighbors */
    if (domain->p.smooth_refine)
    {
        int ng, nb, np;
        int face, corner;
        int nprocs[P4EST_HALF], nblockno, npatchno[P4EST_HALF], nfc;
        int level, max_tlevel;
        int exists;
        int k;
        fclaw2d_patch_relation_t nrel;
        fclaw2d_block_t *block;
        fclaw_patch_t *gpatch, *patch, *npatch;

        /* exchange target refinement level with parallel neighbors */
        p4est_ghost_exchange_custom (wrap->p4est, p4est_wrap_get_ghost (wrap),
                                     sizeof (int),
                                     domain->mirror_target_levels,
                                     domain->ghost_target_levels);

        /* loop through remote patches to process their target level */
        for (ng = 0; ng < domain->num_ghost_patches; ++ng)
        {
            gpatch = domain->ghost_patches + ng;
            gpatch->target_level = domain->ghost_target_levels[ng];
        }

        /* loop through local patches to determine their target level */
        for (nb = 0; nb < domain->num_blocks; ++nb)
        {
            block = domain->blocks + nb;
            for (np = 0; np < block->num_patches; ++np)
            {
                /* compute maximum target level between patch and neighbors */
                patch = block->patches + np;
                level = patch->level;
                max_tlevel = patch->target_level;

                /* loop through face neighbors of this patch */
                for (face = 0; max_tlevel <= level &&
                     face < wrap->p4est_faces; ++face)
                {
                    nrel = fclaw2d_patch_face_neighbors (domain, nb, np, face,
                                                         nprocs, &nblockno,
                                                         npatchno, &nfc);

                    /* we refine ourself if the neighbor wants to be finer */
                    if (nrel == FCLAW2D_PATCH_SAMESIZE)
                    {
                        npatch = fclaw2d_domain_get_neighbor_patch (domain,
                                                                    nprocs[0],
                                                                    nblockno,
                                                                    npatchno
                                                                    [0]);
                        P4EST_ASSERT (npatch->level == level);
                        if (npatch->level >= domain->p.smooth_level)
                            /* Match target level only if we are in a level that
                               should be refined */
                            max_tlevel =
                                SC_MAX (max_tlevel, npatch->target_level);
                    }
                    else if (nrel == FCLAW2D_PATCH_DOUBLESIZE)
                    {
                        npatch = fclaw2d_domain_get_neighbor_patch (domain,
                                                                    nprocs[0],
                                                                    nblockno,
                                                                    npatchno
                                                                    [0]);
                        P4EST_ASSERT (npatch->level == level - 1);
                        if (npatch->level >= domain->p.smooth_level)
                            max_tlevel =
                                SC_MAX (max_tlevel, npatch->target_level);
                    }
                    else if (nrel == FCLAW2D_PATCH_HALFSIZE)
                    {
                        for (k = 0; k < P4EST_HALF; ++k)
                        {
                            npatch =
                                fclaw2d_domain_get_neighbor_patch (domain,
                                                                   nprocs[k],
                                                                   nblockno,
                                                                   npatchno
                                                                   [k]);
                            P4EST_ASSERT (npatch->level == level + 1);
                            if (npatch->level >= domain->p.smooth_level)
                                max_tlevel =
                                    SC_MAX (max_tlevel, npatch->target_level);
                        }
                    }
                    else
                    {
                        FCLAW_ASSERT (nrel == FCLAW2D_PATCH_BOUNDARY);
                    }
                }

#ifdef P4_TO_P8
                /* loop through edge neighbors of this patch */
                /* to be implemented */
#endif

                /* loop through corner neighbors of this patch */
                for (corner = 0; max_tlevel <= level &&
                     corner < wrap->p4est_children; ++corner)
                {
                    exists = fclaw2d_patch_corner_neighbors (domain, nb, np,
                                                             corner, nprocs,
                                                             &nblockno,
                                                             npatchno,
                                                             &nfc, &nrel);
                    if (exists)
                    {
                        FCLAW_ASSERT (nrel != FCLAW2D_PATCH_BOUNDARY);
                        npatch =
                            fclaw2d_domain_get_neighbor_patch (domain,
                                                               nprocs[0],
                                                               nblockno,
                                                               npatchno[0]);
                        if (npatch->level >= domain->p.smooth_level)
                            max_tlevel =
                                SC_MAX (max_tlevel, npatch->target_level);
                    }
                }

                /* pass refinement marker into the plumbing */
                if (max_tlevel < level)
                {
                    FCLAW_ASSERT (max_tlevel == level - 1);
                    p4est_wrap_mark_coarsen (wrap,
                                             (p4est_locidx_t) nb,
                                             (p4est_locidx_t) np);
                }
                else if (max_tlevel > level)
                {
                    /* max target level may be bigger by one or two */
                    FCLAW_ASSERT (max_tlevel <= level + 2);
                    p4est_wrap_mark_refine (wrap,
                                            (p4est_locidx_t) nb,
                                            (p4est_locidx_t) np);
                }
            }
        }
    }

    /* do the adaptation */
    if (p4est_wrap_adapt (wrap))
    {
        fclaw2d_domain_t *newd;

        domain->pp_owned = 0;
        newd = fclaw2d_domain_new (wrap, domain->attributes);
        newd->just_adapted = 1;

        fclaw2d_domain_copy_parameters (newd, domain);
        return newd;
    }
    else
    {
        /* clean up the target levels since we don't adapt */
        int nb, np;
        fclaw2d_block_t *block;
        fclaw_patch_t *patch;

        for (nb = 0; nb < domain->num_blocks; ++nb)
        {
            block = domain->blocks + nb;
            for (np = 0; np < block->num_patches; ++np)
            {
                /* compute maximum target level between patch and neighbors */
                patch = block->patches + np;
                patch->target_level = patch->level;
            }
        }
        return NULL;
    }
}

fclaw2d_domain_t *
fclaw2d_domain_partition (fclaw2d_domain_t * domain, int weight_exponent)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_locidx_t uf, ul, uof;

    FCLAW_ASSERT (domain->pp_owned);
    FCLAW_ASSERT (domain->just_adapted);
    FCLAW_ASSERT (!domain->just_partitioned);

    domain->just_adapted = 0;
    if (p4est_wrap_partition (wrap, weight_exponent, &uf, &ul, &uof))
    {
        fclaw2d_domain_t *newd;

        domain->pp_owned = 0;
        newd = fclaw2d_domain_new (wrap, domain->attributes);
        newd->just_partitioned = 1;
        newd->partition_unchanged_first = (int) uf;
        newd->partition_unchanged_length = (int) ul;
        newd->partition_unchanged_old_first = (int) uof;

        fclaw2d_domain_copy_parameters (newd, domain);
        return newd;
    }
    else
    {
        return NULL;
    }
}

void
fclaw2d_domain_partition_unchanged (fclaw2d_domain_t * domain,
                                    int *unchanged_first,
                                    int *unchanged_length,
                                    int *unchanged_old_first)
{
    FCLAW_ASSERT (domain->pp_owned);
    FCLAW_ASSERT (!domain->just_adapted);
    FCLAW_ASSERT (domain->just_partitioned);

    if (unchanged_first != NULL)
    {
        *unchanged_first = domain->partition_unchanged_first;
    }
    if (unchanged_length != NULL)
    {
        *unchanged_length = domain->partition_unchanged_length;
    }
    if (unchanged_old_first != NULL)
    {
        *unchanged_old_first = domain->partition_unchanged_old_first;
    }
}

void
fclaw2d_domain_complete (fclaw2d_domain_t * domain)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (domain->pp_owned);
    FCLAW_ASSERT (!domain->just_adapted);
    FCLAW_ASSERT (domain->just_partitioned);

    domain->just_partitioned = 0;
    domain->partition_unchanged_first = 0;
    domain->partition_unchanged_length = 0;
    domain->partition_unchanged_old_first = 0;

    p4est_wrap_complete (wrap);
}

void
fclaw2d_domain_write_vtk (fclaw2d_domain_t * domain, const char *basename)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (wrap != NULL);
    FCLAW_ASSERT (wrap->p4est != NULL);

    p4est_vtk_write_file (wrap->p4est, NULL, basename);
}

static void
fclaw2d_domain_list_level_callback (fclaw2d_domain_t * domain,
                                    fclaw_patch_t * patch, int block_no,
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
                                        fclaw_patch_t * patch, int block_no,
                                        int patch_no, void *user)
{
    fclaw2d_domain_list_neighbors_t *ln =
        (fclaw2d_domain_list_neighbors_t *) user;
    fclaw2d_patch_relation_t fnt;
    int faceno, cornerno, rcorner;
    int rproc[P4EST_FACES], rblockno, rpatchno[P4EST_FACES], rfaceno;

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
#ifdef P4_TO_P8
    /* to be implemented: list edge neighbors as well */
#endif
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
                                      fclaw_patch_t * old_patch,
                                      fclaw2d_domain_t * new_domain,
                                      fclaw_patch_t * new_patch,
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

typedef struct fclaw2d_search_data
{
    fclaw2d_domain_t *domain;
    sc_array_t *coordinates;
    sc_array_t *results;
}
fclaw2d_search_data_t;

static int
search_point_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                 void *point)
{
    int ip, jb;
    int earlier, now;
    int *pentry;
    double x, y;
#ifdef P4_TO_P8
    double z;
#endif
    double *xyentry;
    fclaw2d_block_t *block;
#ifdef FCLAW_ENABLE_DEBUG
    fclaw_patch_t *patch;
#endif
    fclaw2d_search_data_t *sd = (fclaw2d_search_data_t *) p4est->user_pointer;
    p4est_qcoord_t qh;

    FCLAW_ASSERT (sd != NULL);
    FCLAW_ASSERT (sd->domain != NULL);
    FCLAW_ASSERT (sd->coordinates != NULL);
    FCLAW_ASSERT (sd->results != NULL);
    FCLAW_ASSERT (point != NULL);

    /* access point data */
    pentry = (int *) point;
    jb = pentry[0];
    if (jb != (int) which_tree)
    {
        /* this is the wrong tree entirely */
        return 0;
    }
    ip = pentry[1];
    xyentry = (double *) sc_array_index_int (sd->coordinates, ip);

    /* compare quadrant coordinates with point coordinates */
    qh = P4EST_QUADRANT_LEN (quadrant->level);
    x = xyentry[0];
    if (x < quadrant->x * fclaw2d_smallest_h ||
        x > (quadrant->x + qh) * fclaw2d_smallest_h)
    {
        return 0;
    }
    y = xyentry[1];
    if (y < quadrant->y * fclaw2d_smallest_h ||
        y > (quadrant->y + qh) * fclaw2d_smallest_h)
    {
        return 0;
    }
#ifdef P4_TO_P8
    z = xyentry[2];
    if (z < quadrant->z * fclaw2d_smallest_h ||
        z > (quadrant->z + qh) * fclaw2d_smallest_h)
    {
        return 0;
    }
#endif
    /* now we know that the point is contained in the search quadrant */

    if (local_num < 0)
    {
        FCLAW_ASSERT (local_num == -1);
        /* this is not a leaf */
        return 1;
    }
    /* now we know that we have found the point in a leaf */

    FCLAW_ASSERT (0 <= jb && jb < sd->domain->num_blocks);
    block = sd->domain->blocks + jb;
    FCLAW_ASSERT (block != NULL);
    now = (int) local_num - block->num_patches_before;
    FCLAW_ASSERT (0 <= now && now < block->num_patches);
#ifdef FCLAW_ENABLE_DEBUG
    patch = block->patches + now;
    FCLAW_ASSERT (patch != NULL);
#endif

    /* do the check a second time with the patch data */
    FCLAW_ASSERT (x >= patch->d2->xlower && x <= patch->d2->xupper);
    FCLAW_ASSERT (y >= patch->d2->ylower && y <= patch->d2->yupper);
#ifdef P4_TO_P8
    FCLAW_ASSERT (z >= patch->d3->zlower && z <= patch->d3->zupper);
#endif

    /* remember the smallest local quadrant number as result */
    earlier = *(int *) sc_array_index_int (sd->results, ip);
    FCLAW_ASSERT (-1 <= earlier && earlier < block->num_patches);
    if (earlier >= 0)
    {
        now = SC_MIN (earlier, now);
    }
    if (now != earlier)
    {
        P4EST_ASSERT (earlier == -1 || now < earlier);
        *(int *) sc_array_index_int (sd->results, ip) = now;
    }

    /* for leaves the return value is irrelevant */
    return 1;
}

void
fclaw2d_domain_search_points (fclaw2d_domain_t * domain,
                              sc_array_t * block_offsets,
                              sc_array_t * coordinates, sc_array_t * results)
{
    int ip, jb, mpiret;
    int num_blocks;
    int num_points;
    int pbegin, pend;
    int *pentry;
    int *found;
    int *found_buffer;
    int *resi;
    sc_array_t *points;

    p4est_t *p4est;
    p4est_wrap_t *wrap;
    void *user_save;

    fclaw2d_search_data_t search_data, *sd = &search_data;

    /* assert validity of parameters */
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (block_offsets != NULL);
    FCLAW_ASSERT (coordinates != NULL);
    FCLAW_ASSERT (results != NULL);

    num_blocks = domain->num_blocks;
    FCLAW_ASSERT (0 <= num_blocks);

    FCLAW_ASSERT (block_offsets->elem_size == sizeof (int));
    FCLAW_ASSERT (coordinates->elem_size == P4EST_DIM * sizeof (double));
    FCLAW_ASSERT (results->elem_size == sizeof (int));

    FCLAW_ASSERT (num_blocks + 1 == (int) block_offsets->elem_count);
    FCLAW_ASSERT (0 == *(int *) sc_array_index_int (block_offsets, 0));
    num_points = *(int *) sc_array_index_int (block_offsets, num_blocks);

    FCLAW_ASSERT (num_points == (int) coordinates->elem_count);
    FCLAW_ASSERT (num_points == (int) results->elem_count);

    /* prepare results array */
    for (ip = 0; ip < num_points; ++ip)
    {
        *(int *) sc_array_index_int (results, ip) = -1;
    }

    /* construct input set for p4est search: block and point number */
    points = sc_array_new_size (2 * sizeof (int), num_points);
    pbegin = 0;
    for (jb = 0; jb < num_blocks; ++jb)
    {
        pend = *(int *) sc_array_index_int (block_offsets, jb + 1);
        FCLAW_ASSERT (pbegin <= pend);

        /* managemant data used internal to the p4est search */
        for (ip = pbegin; ip < pend; ++ip)
        {
            pentry = (int *) sc_array_index (points, ip);
            pentry[0] = jb;
            pentry[1] = ip;
        }
        pbegin = pend;
    }
    FCLAW_ASSERT (pbegin == num_points);

    /* stash relevant information to pass to search */
    sd->domain = domain;
    sd->coordinates = coordinates;
    sd->results = results;

    /* process-local search through p4est */
    wrap = (p4est_wrap_t *) domain->pp;
    FCLAW_ASSERT (wrap != NULL);
    p4est = wrap->p4est;
    FCLAW_ASSERT (p4est != NULL);
    FCLAW_ASSERT (p4est->connectivity != NULL);
    FCLAW_ASSERT (p4est->connectivity->num_trees ==
                  (p4est_topidx_t) num_blocks);
    user_save = p4est->user_pointer;
    p4est->user_pointer = sd;
    p4est_search_local (p4est, 0, NULL, search_point_fn, points);
    p4est->user_pointer = user_save;

    /* synchronize results in parallel */
    found = FCLAW_ALLOC (int, num_points);
    for (ip = 0; ip < num_points; ++ip)
    {
        resi = (int *) sc_array_index (results, ip);
        found[ip] = (*resi < 0 ? domain->mpisize : domain->mpirank);
    }

    found_buffer = FCLAW_ALLOC (int, num_points);
    mpiret = sc_MPI_Allreduce (found, found_buffer, num_points, sc_MPI_INT,
                               sc_MPI_MIN, domain->mpicomm);
    SC_CHECK_MPI (mpiret);
    FCLAW_FREE (found);

    for (ip = 0; ip < num_points; ++ip)
    {
        if (found_buffer[ip] != domain->mpirank) {
            *(int *) sc_array_index (results, ip) = -1;
        }
    }

    /* tidy up memory */
    FCLAW_FREE (found_buffer);
    sc_array_destroy (points);
}

typedef struct fclaw2d_ray_integral
{
    void *ray;
    double *integral;
}
fclaw2d_ray_integral_t;

typedef struct fclaw2d_integrate_ray_data
{
    fclaw2d_domain_t *domain;
    fclaw2d_integrate_ray_t integrate_ray;
    void *user;
}
fclaw2d_integrate_ray_data_t;

static int
integrate_ray_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                  void *point)
{
    double integral;
    int intersects;
    fclaw2d_domain_t *domain;
    fclaw_patch_t *patch;


    fclaw_patch_t fclaw2d_patch;
    memset (&fclaw2d_patch, 0, sizeof (fclaw_patch_t));

    fclaw_patch_bounds_2d_t fclaw2d_patch_bounds;
    fclaw2d_patch.d2 = &fclaw2d_patch_bounds;

    int patchno;

    /* assert that the user_pointer contains a valid integrate_ray_data_t */
    fclaw2d_integrate_ray_data_t *ird
        = (fclaw2d_integrate_ray_data_t *) p4est->user_pointer;
    FCLAW_ASSERT (ird != NULL);
    FCLAW_ASSERT (ird->domain != NULL);
    FCLAW_ASSERT (ird->integrate_ray != NULL);

    /* assert that point is a valid ray_integral_t */
    fclaw2d_ray_integral_t *ri = (fclaw2d_ray_integral_t *) point;
    FCLAW_ASSERT (ri != NULL);
    FCLAW_ASSERT (ri->ray != NULL);
    FCLAW_ASSERT (ri->integral != NULL);

    /* collect patch information */
    domain = ird->domain;
    if (local_num >= 0)
    {
        fclaw2d_block_t *block = domain->blocks + which_tree;
        patchno = local_num - block->num_patches_before;
        patch = block->patches + patchno;
    }
    else
    {
        /* create artifical patch and fill it based on the quadrant */
        patchno = -1;
        patch = &fclaw2d_patch;
        patch->level = quadrant->level;
        patch->target_level = quadrant->level;
        patch->flags = p4est_quadrant_child_id (quadrant);
        fclaw2d_patch_set_boundary_xylower (patch, quadrant);
        patch->u.next = NULL;
        patch->user = NULL;
    }

    /* compute local integral and add it onto the ray integral */
    integral = 0.;
    intersects = ird->integrate_ray (domain, patch, which_tree, patchno,
                                     ri->ray, &integral, ird->user);
    if (local_num >= 0)
    {
        *(ri->integral) += integral;
    }
    return intersects;
}

void
fclaw2d_domain_integrate_rays (fclaw2d_domain_t * domain,
                               fclaw2d_integrate_ray_t intersect,
                               sc_array_t * rays, sc_array_t * integrals,
                               void * user)
{
    int i;
    size_t nintz;
    sc_array_t lints[1];
    sc_array_t ri[1];
    fclaw2d_ray_integral_t *rayint;
    fclaw2d_integrate_ray_data_t integrate_ray_data, *ird =
        &integrate_ray_data;
    p4est_t p4est;
    p4est_wrap_t *wrap;

    /* assert validity of parameters */
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (intersect != NULL);
    FCLAW_ASSERT (rays != NULL);
    FCLAW_ASSERT (integrals != NULL);
    FCLAW_ASSERT (integrals->elem_size == sizeof (double));
    FCLAW_ASSERT (rays->elem_count == integrals->elem_count);

    /* create local storage for integral values */
    nintz = integrals->elem_count;
    sc_array_init_count (lints, sizeof (double), nintz);
    memset (lints->array, 0, sizeof (double) * nintz);

    /* construct ray_integral_t array from rays */
    sc_array_init_count (ri, sizeof (fclaw2d_ray_integral_t), nintz);
    for (i = 0; i < (int) nintz; ++i)
    {
        rayint = (fclaw2d_ray_integral_t *) sc_array_index_int (ri, i);
        rayint->ray = sc_array_index_int (rays, i);
        rayint->integral = (double *) sc_array_index_int (lints, i);
    }

    /* construct fclaw2d_integrate_ray_data_t */
    ird->domain = domain;
    ird->integrate_ray = intersect;
    ird->user = user;

    /* process-local integration through p4est */
    wrap = (p4est_wrap_t *) domain->pp;
    FCLAW_ASSERT (wrap != NULL);
    FCLAW_ASSERT (wrap->p4est != NULL);
    p4est = *(wrap->p4est);
    p4est.user_pointer = ird;
    p4est_search_local (&p4est, 0, NULL, integrate_ray_fn, ri);

    /* allreduce local integral values in parallel */
    sc_MPI_Allreduce (lints->array, integrals->array, nintz,
                      sc_MPI_DOUBLE, sc_MPI_SUM, domain->mpicomm);

    sc_array_reset (ri);
    sc_array_reset (lints);
}

/******* code for overlap exchange algorithm *******/

typedef enum comm_tag
{
    COMM_TAG_CONSDATA = 5526,
    COMM_TAG_PRODATA = 5527
}
comm_tag_t;

typedef struct overlap_query_ind
{
    int rank;
    sc_array_t oqs;
}
overlap_ind_t;

typedef struct overlap_send_buf
{
    int rank;
    sc_array_t ops;
}
overlap_buf_t;

typedef struct overlap_point
{
    void *point;
    int rank;
    size_t id;
}
overlap_point_t;

typedef struct overlap_producer_comm
{
    fclaw2d_domain_t *domain;
    fclaw2d_interpolate_point_t interpolate;
    void *user;
    p4est_t *pro4est;
    sc_MPI_Comm glocomm;
    int prorank;
    int iprorank;
    size_t point_size;
    sc_array_t *recv_buffer;
    sc_array_t *recv_reqs;
    sc_array_t *send_reqs;
}
overlap_producer_comm_t;

typedef struct overlap_consumer_comm
{
    fclaw2d_domain_t *domain;
    fclaw2d_interpolate_point_t interpolate;
    void *user;
    sc_array_t *query_points;
    sc_array_t *query_indices;
    sc_MPI_Comm glocomm;
    int conrank;
    int iconrank;
    size_t point_size;
    sc_array_t *send_buffer;
    sc_array_t *send_reqs;
    sc_array_t *recv_buffer;
    sc_array_t *recv_reqs;
}
overlap_consumer_comm_t;

typedef struct overlap_global_comm
{
    sc_MPI_Comm glocomm;
    overlap_producer_comm_t pro, *p;
    overlap_consumer_comm_t con, *c;
}
overlap_global_comm_t;

static void
overlap_consumer_add (overlap_consumer_comm_t * c, void *point, int rank)
{
    size_t bcount;
    overlap_buf_t *sb;
    overlap_ind_t *qi;

    P4EST_ASSERT (c != NULL);
    P4EST_ASSERT (c->send_buffer != NULL && c->query_indices != NULL);
    P4EST_ASSERT (0 <= rank && rank < c->domain->mpisize);
    overlap_point_t *op = (overlap_point_t *) point;
    FCLAW_ASSERT (op != NULL && op->point != NULL);
    op->rank = rank;            /* mark, that we added this point to the process buffer */

    /* if we have a new rank, push new send buffer */
    bcount = c->send_buffer->elem_count;
    sb = NULL;
    qi = NULL;
    if (bcount > 0)
    {
        sb = (overlap_buf_t *) sc_array_index (c->send_buffer, bcount - 1);
        qi = (overlap_ind_t *) sc_array_index (c->query_indices, bcount - 1);
        P4EST_ASSERT (sb->rank == qi->rank);
        P4EST_ASSERT (sb->rank <= rank);
        P4EST_ASSERT (sb->ops.elem_count == qi->oqs.elem_count);
        P4EST_ASSERT (sb->ops.elem_count > 0);
    }
    if (bcount == 0 || sb->rank < rank)
    {
        sb = (overlap_buf_t *) sc_array_push (c->send_buffer);
        qi = (overlap_ind_t *) sc_array_push (c->query_indices);
        sb->rank = qi->rank = rank;
        sc_array_init (&sb->ops, c->point_size);
        sc_array_init (&qi->oqs, sizeof (size_t));
    }
    memcpy (sc_array_push (&sb->ops), op->point, c->point_size);
    memcpy (sc_array_push (&qi->oqs), &op->id, qi->oqs.elem_size);
}

static int
interpolate_partition_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t * quadrant, int pfirst, int plast,
                          void *point)
{
    fclaw2d_domain_t *domain;
    fclaw_patch_t *patch;

    fclaw_patch_t fclaw2d_patch;
    memset(&fclaw2d_patch, 0, sizeof(fclaw_patch_t));

    fclaw_patch_bounds_2d_t fclaw2d_patch_bounds;
    fclaw2d_patch.d2 = &fclaw2d_patch_bounds;

    int patchno;
    int intersects;

    /* assert that the point is a valid overlap_point_t and was not added yet */
    overlap_point_t *op = (overlap_point_t *) point;
    FCLAW_ASSERT (op != NULL && op->point != NULL);
    if (op->rank >= 0)
    {
        return 0;               /* avoid sending the same point twice */
    }

    /* assert that the user_pointer contains a valid interpolation_data_t */
    overlap_consumer_comm_t *c
        = (overlap_consumer_comm_t *) p4est->user_pointer;
    FCLAW_ASSERT (c != NULL);
    FCLAW_ASSERT (c->domain != NULL && c->domain->pp != NULL);
    FCLAW_ASSERT (c->interpolate != NULL);

    /* create artifical domain, that only contains mpi and tree structure data */
    domain = FCLAW_ALLOC (fclaw2d_domain_t, 1);
    /* Todo: works only for congruent communicators. Do we really need all this
     * information? */
    fclaw2d_domain_init_meta (domain, (pfirst == plast) ? pfirst : -1);

    /* create artifical patch and fill it based on the quadrant */
    patch = &fclaw2d_patch;
    patch->level = quadrant->level;
    patch->target_level = quadrant->level;
    patch->flags = p4est_quadrant_child_id (quadrant);
    fclaw2d_patch_set_boundary_xylower (patch, quadrant);
    patch->u.next = NULL;
    patch->user = NULL;
    patchno = -1;               /* marks patch as artifical patch */

    /* check if the patch (or its decendants) may contribute to the
     * interpolation data of the point */
    intersects =
        c->interpolate (domain, patch, which_tree, patchno, op->point,
                        c->user);

    fclaw2d_domain_destroy (domain);

    if (!intersects)
    {
        return 0;
    }

    /* we have located the point in the intersection quadrant */
    if (pfirst == plast)
    {
        /* we have intersected with a leaf quadrant */
        overlap_consumer_add (c, op, pfirst);
    }
    return 1;
}

static int
interpolate_local_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                      p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                      void *point)
{
    fclaw2d_domain_t *domain;
    fclaw_patch_t *patch;

    fclaw_patch_t fclaw2d_patch;
    memset(&fclaw2d_patch, 0, sizeof(fclaw_patch_t));

    fclaw_patch_bounds_2d_t fclaw2d_patch_bounds;
    fclaw2d_patch.d2 = &fclaw2d_patch_bounds;

    int patchno;

    /* assert that the user_pointer contains a valid interpolation_data_t */
    overlap_producer_comm_t *p
        = (overlap_producer_comm_t *) p4est->user_pointer;
    FCLAW_ASSERT (p != NULL);
    FCLAW_ASSERT (p->domain != NULL);
    FCLAW_ASSERT (p->interpolate != NULL);

    domain = p->domain;
    if (local_num >= 0)
    {
        fclaw2d_block_t *block = domain->blocks + which_tree;
        patchno = local_num - block->num_patches_before;
        patch = block->patches + patchno;
    }
    else
    {
        /* create artifical patch and fill it based on the quadrant */
        patchno = -1;
        patch = &fclaw2d_patch;
        patch->level = quadrant->level;
        patch->target_level = quadrant->level;
        patch->flags = p4est_quadrant_child_id (quadrant);
        fclaw2d_patch_set_boundary_xylower (patch, quadrant);
        patch->u.next = NULL;
        patch->user = NULL;
    }

    return p->interpolate (domain, patch, which_tree, patchno, point,
                           p->user);
}

#ifdef FCLAW_ENABLE_MPI
static void
consumer_producer_notify (overlap_global_comm_t * g)
{
    overlap_producer_comm_t *p = g->p;
    overlap_consumer_comm_t *c = g->c;
    size_t bz, bcount;
    sc_array_t *receivers, *senders;
    sc_array_t *payload_in, *payload_out;
    int num_receivers, num_senders;
    overlap_buf_t *sb, *rb;
    int same_rank, num_ops, i;
    int mpiret;

    /* assemble and execute receiver and payload query */
    num_receivers = (int) (bcount = c->send_buffer->elem_count);
    receivers = sc_array_new_count (sizeof (int), bcount);
    senders = sc_array_new (sizeof (p4est_locidx_t));
    payload_in = sc_array_new_count (sizeof (int), bcount);
    payload_out = sc_array_new (sizeof (p4est_locidx_t));
    for (bz = 0; bz < bcount; ++bz)
    {
        sb = (overlap_buf_t *) sc_array_index (c->send_buffer, bz);
        *(int *) sc_array_index (receivers, bz) = sb->rank;
        *(p4est_locidx_t *) sc_array_index (payload_in, bz) =
            (p4est_locidx_t) sb->ops.elem_count;
    }
    sc_notify_ext (receivers, senders, payload_in, payload_out, g->glocomm);
    num_senders = (int) senders->elem_count;
    P4EST_LDEBUGF ("Overlap exchange receivers %d senders %d\n",
                   num_receivers, num_senders);

    /* post nonblocking receives for the point data of the consumer side */
    p->recv_buffer = sc_array_new_count (sizeof (overlap_buf_t), num_senders);
    p->recv_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
    p->iprorank = c->iconrank = -1;
    for (i = 0; i < num_senders; ++i)
    {
        /* initalize and allocate the buffer according to the payload */
        rb = (overlap_buf_t *) sc_array_index_int (p->recv_buffer, i);
        rb->rank = *(int *) sc_array_index_int (senders, i);
        same_rank = (rb->rank == p->prorank);
        num_ops =
            same_rank ? 0 : *(int *) sc_array_index_int (payload_out, i);
        sc_array_init_size (&(rb->ops), p->point_size, (size_t) num_ops);
        if (same_rank)
        {
            p->iprorank = i;    /* save the index in the producer buffer */
            *(sc_MPI_Request *) sc_array_index_int (p->recv_reqs, i) =
                sc_MPI_REQUEST_NULL;
            continue;
        }

        /* receive the array of overlap_point_t data and store it in the buffer */
        mpiret =
            sc_MPI_Irecv (rb->ops.array, num_ops * p->point_size,
                          sc_MPI_BYTE, rb->rank, COMM_TAG_CONSDATA,
                          p->glocomm,
                          (sc_MPI_Request *) sc_array_index_int (p->recv_reqs,
                                                                 i));
        SC_CHECK_MPI (mpiret);
    }

    sc_array_destroy (receivers);
    sc_array_destroy (senders);
    sc_array_destroy (payload_in);
    sc_array_destroy (payload_out);
}

static void
consumer_post_messages (overlap_consumer_comm_t * c)
{
    overlap_buf_t *sb, *rb;
    int num_receivers, same_rank, num_ops, i;
    int mpiret;

    /* send the point data to the producer side in a nonblocking way */
    num_receivers = (int) c->send_buffer->elem_count;
    c->send_reqs =
        sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i)
    {
        sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, i);

        if (sb->rank == c->conrank)
        {
            c->iconrank = i;    /* save the index in the consumer buffer */
            *(sc_MPI_Request *) sc_array_index_int (c->send_reqs, i) =
                sc_MPI_REQUEST_NULL;
            continue;
        }

        mpiret =
            sc_MPI_Isend (sb->ops.array,
                          sb->ops.elem_count * c->point_size,
                          sc_MPI_BYTE, sb->rank, COMM_TAG_CONSDATA,
                          c->glocomm,
                          (sc_MPI_Request *) sc_array_index_int (c->send_reqs,
                                                                 i));
        SC_CHECK_MPI (mpiret);
    }

    /* recv the updated point data from the producer side in a nonblocking way */
    c->recv_reqs =
        sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    c->recv_buffer =
        sc_array_new_size (sizeof (overlap_buf_t), num_receivers);
    for (i = 0; i < num_receivers; ++i)
    {
        rb = (overlap_buf_t *) sc_array_index_int (c->recv_buffer, i);
        sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, i);
        rb->rank = sb->rank;
        same_rank = (rb->rank == c->conrank);
        num_ops = same_rank ? 0 : (int) sb->ops.elem_count;
        sc_array_init_size (&(rb->ops), c->point_size, (size_t) num_ops);

        if (same_rank)
        {
            *(sc_MPI_Request *) sc_array_index_int (c->recv_reqs,
                                                    c->iconrank) =
                sc_MPI_REQUEST_NULL;
            continue;
        }

        /* receive the array of overlap_point_t data and store it in the buffer */
        mpiret =
            sc_MPI_Irecv (rb->ops.array,
                          rb->ops.elem_count * c->point_size,
                          sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, c->glocomm,
                          (sc_MPI_Request *) sc_array_index_int (c->recv_reqs,
                                                                 i));
        SC_CHECK_MPI (mpiret);
    }
}

static void
producer_interpolate (overlap_producer_comm_t * p)
{
    overlap_buf_t *rb;
    int num_senders, i;
    int remaining, received;
    int *prod_indices;
    int mpiret;

    /* compute producer data for all incoming messages as soon as they come in */
    num_senders = (int) p->recv_reqs->elem_count;
    prod_indices = FCLAW_ALLOC (int, num_senders);
    p->send_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
    remaining = num_senders;
    if (p->iprorank >= 0)
    {
        *(sc_MPI_Request *) sc_array_index_int (p->send_reqs, p->iprorank) =
            sc_MPI_REQUEST_NULL;
        remaining--;            /* since we set the iprorank-th request to null earlier */
    }
    while (remaining > 0)
    {
        mpiret =
            sc_MPI_Waitsome (num_senders,
                             (sc_MPI_Request *) p->recv_reqs->array,
                             &received, prod_indices, sc_MPI_STATUSES_IGNORE);
        SC_CHECK_MPI (mpiret);
        P4EST_ASSERT (received != sc_MPI_UNDEFINED);
        P4EST_ASSERT (received > 0);

        for (i = 0; i < received; ++i)
        {
            /* compute the prodata for the points sent by process prod_indices[i] */
            P4EST_ASSERT (0 <= prod_indices[i]
                          && prod_indices[i] < num_senders);
            rb = (overlap_buf_t *) sc_array_index_int (p->recv_buffer,
                                                       prod_indices[i]);
            p4est_search_local (p->pro4est, 0, NULL, interpolate_local_fn,
                                &(rb->ops));

            /* send the requested producer data back in a nonblocking way */
            mpiret =
                sc_MPI_Isend (rb->ops.array,
                              rb->ops.elem_count * p->point_size,
                              sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA,
                              p->glocomm,
                              (sc_MPI_Request *)
                              sc_array_index_int (p->send_reqs,
                                                  prod_indices[i]));
            SC_CHECK_MPI (mpiret);
        }

        remaining -= received;
    }

    FCLAW_FREE (prod_indices);
}
#endif /* FCLAW_ENABLE_MPI */

static void
consumer_update_from_buffer (overlap_consumer_comm_t * c, sc_array_t * buffer,
                             int bi)
{
    overlap_buf_t *rb;
    overlap_ind_t *qi;
    void *p, *qp;
    size_t *pi;
    int i;

    /* obtain the array of points we want to update the query points with */
    P4EST_ASSERT (0 <= bi && bi < (int) buffer->elem_count);
    rb = (overlap_buf_t *) sc_array_index_int (buffer, bi);
    qi = (overlap_ind_t *) sc_array_index_int (c->query_indices, bi);

    /* copy updated points into the query-point array */
    for (i = 0; i < (int) rb->ops.elem_count; ++i)
    {
        p = sc_array_index_int (&rb->ops, i);
        pi = (size_t *) sc_array_index_int (&qi->oqs, i);
        qp = sc_array_index (c->query_points, *pi);
        memcpy (qp, p, rb->ops.elem_size);
    }
}

#ifdef FCLAW_ENABLE_MPI
static void
consumer_update_query_points (overlap_consumer_comm_t * c)
{
    int num_receivers, i;
    int remaining, received;
    int *cons_indices;
    int mpiret;

    /* compute producer data for all incoming messages as soon as they come in */
    num_receivers = (int) c->recv_reqs->elem_count;
    cons_indices = FCLAW_ALLOC (int, num_receivers);
    remaining = (c->iconrank >= 0) ? num_receivers - 1 : num_receivers;
    while (remaining > 0)
    {
        mpiret =
            sc_MPI_Waitsome (num_receivers,
                             (sc_MPI_Request *) c->recv_reqs->array,
                             &received, cons_indices, sc_MPI_STATUSES_IGNORE);
        SC_CHECK_MPI (mpiret);
        P4EST_ASSERT (received != sc_MPI_UNDEFINED);
        P4EST_ASSERT (received > 0);

        for (i = 0; i < received; ++i)
        {
            consumer_update_from_buffer (c, c->recv_buffer, cons_indices[i]);
        }

        remaining -= received;
    }

    FCLAW_FREE (cons_indices);
}

static void
consumer_waitall (overlap_consumer_comm_t * c)
{
    int mpiret;
    int num_receivers;

    /* wait for the nonblocking sends to complete */
    num_receivers = (int) c->send_reqs->elem_count;
    mpiret =
        sc_MPI_Waitall (num_receivers, (sc_MPI_Request *) c->send_reqs->array,
                        sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
}

static void
producer_waitall (overlap_producer_comm_t * p)
{
    int mpiret;
    int num_senders;

    /* wait for the nonblocking sends to complete */
    num_senders = (int) p->send_reqs->elem_count;
    mpiret =
        sc_MPI_Waitall (num_senders, (sc_MPI_Request *) p->send_reqs->array,
                        sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
}
#endif

static void
consumer_free_communication_data (overlap_consumer_comm_t * c)
{
    overlap_ind_t *qi;
    overlap_buf_t *sb;
    int i, num_queries;
#ifdef FCLAW_ENABLE_MPI
    overlap_buf_t *rb;
#ifdef FCLAW_ENABLE_DEBUG
    int prev_rank;
#endif
    size_t bz, bcount;

    sc_array_destroy (c->recv_reqs);
    sc_array_destroy (c->send_reqs);
#ifdef FCLAW_ENABLE_DEBUG
    prev_rank = -1;
#endif
    bcount = c->send_buffer->elem_count;
    for (bz = 0; bz < bcount; ++bz)
    {
        sb = (overlap_buf_t *) sc_array_index (c->send_buffer, bz);
        rb = (overlap_buf_t *) sc_array_index (c->recv_buffer, bz);
        SC_ASSERT (sb->rank == rb->rank);
#ifdef FCLAW_ENABLE_DEBUG
        SC_ASSERT (prev_rank < sb->rank);
        prev_rank = sb->rank;
        if (bz == (size_t) c->iconrank)
        {
            P4EST_ASSERT (rb->ops.elem_count == 0);
        }
        else
        {
            P4EST_ASSERT (rb->ops.elem_count == sb->ops.elem_count);
        }
#endif
        FCLAW_ASSERT (sb->ops.elem_count > 0);
        sc_array_reset (&sb->ops);
        sc_array_reset (&rb->ops);
    }
    sc_array_destroy_null (&c->recv_buffer);
#else /* !FCLAW_ENABLE_MPI */
    if (c->send_buffer->elem_count)
    {
        sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, 0);
        sc_array_reset (&sb->ops);
    }
#endif
    num_queries = (int) c->query_indices->elem_count;
    for (i = 0; i < num_queries; i++)
    {
        qi = (overlap_ind_t *) sc_array_index_int (c->query_indices, i);
        sc_array_reset (&qi->oqs);
    }
    sc_array_destroy_null (&c->query_indices);
    sc_array_destroy_null (&c->send_buffer);
}

static void
producer_free_communication_data (overlap_producer_comm_t * p)
{
#ifdef FCLAW_ENABLE_MPI
    overlap_buf_t *rb;
    int num_senders, i;

    sc_array_destroy (p->recv_reqs);
    sc_array_destroy (p->send_reqs);
    num_senders = (int) p->recv_buffer->elem_count;
    for (i = 0; i < num_senders; ++i)
    {
        rb = (overlap_buf_t *) sc_array_index_int (p->recv_buffer, i);
#ifdef FCLAW_ENABLE_DEBUG
        if (i == p->iprorank)
        {
            P4EST_ASSERT (rb->ops.elem_count == 0);
        }
        else
        {
            P4EST_ASSERT (rb->ops.elem_count > 0);
        }
#endif
        sc_array_reset (&rb->ops);
    }
    sc_array_destroy_null (&p->recv_buffer);
#endif
}

static void
consumer_producer_update_local (overlap_global_comm_t * g)
{
    overlap_consumer_comm_t *c = g->c;
    overlap_producer_comm_t *p = g->p;
    overlap_buf_t *sb;

    if (c->iconrank >= 0 && c->send_buffer->elem_count)
    {
        /* Interpolate point-data of local points. Instead of copying to the
         * producer buffer, we update the points in-place. */
        sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer,
                                                   c->iconrank);
        p4est_search_local (p->pro4est, 0, NULL, interpolate_local_fn,
                            &(sb->ops));
        consumer_update_from_buffer (c, c->send_buffer, c->iconrank);
    }
}

void
fclaw2d_overlap_exchange (fclaw2d_domain_t * domain,
                          sc_array_t * query_points,
                          fclaw2d_interpolate_point_t interpolate, void *user)
{
    p4est_t p4est;
    p4est_connectivity_t conn;
    p4est_wrap_t *wrap;
    size_t iz, nipz;
    sc_array_t *iquery_points;
    overlap_point_t *ip;
    overlap_global_comm_t global, *g = &global;
    overlap_producer_comm_t *p = g->p = &g->pro;
    overlap_consumer_comm_t *c = g->c = &g->con;

    /* assert validity of parameters */
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (query_points != NULL);
    FCLAW_ASSERT (interpolate != NULL);

    fclaw_global_essentialf ("OVERLAP: exchange partition\n");

    /* extract p4est data from wrapper */
    wrap = (p4est_wrap_t *) domain->pp;
    FCLAW_ASSERT (wrap != NULL);
    FCLAW_ASSERT (wrap->p4est != NULL);
    FCLAW_ASSERT (wrap->conn != NULL);
    p4est = *(wrap->p4est);
    conn = *(wrap->conn);

    /* initialize communication data */
    /*we assume that consumer and producer operate on congruent communicators */
    g->glocomm = c->glocomm = p->glocomm = domain->mpicomm;
    c->domain = p->domain = domain;
    c->interpolate = p->interpolate = interpolate;
    c->user = p->user = user;
    c->query_points = query_points;
    c->query_indices = sc_array_new (sizeof (overlap_ind_t));
    c->point_size = p->point_size = query_points->elem_size;
    c->conrank = p->prorank = domain->mpirank;
    c->send_buffer = sc_array_new (sizeof (overlap_buf_t));
    p4est.user_pointer = p;
    p->pro4est = &p4est;

    fclaw_global_essentialf ("OVERLAP: customer partition search\n");

    /* store the query points as overlap_point_t to be able to mark their last
     * appearance in the search */
    nipz = query_points->elem_count;
    iquery_points = sc_array_new_count (sizeof (overlap_point_t), nipz);
    for (iz = 0; iz < nipz; ++iz)
    {
        ip = (overlap_point_t *) sc_array_index (iquery_points, iz);
        ip->point = sc_array_index (query_points, iz);
        ip->rank = -1;
        ip->id = iz;
    }

    /* search for the query points in the producer-partition and create a buffer
     * to send them to the respective producer ranks */
    p4est_search_partition_gfx (p4est.global_first_quadrant,
                                p4est.global_first_position, p4est.mpisize,
                                conn.num_trees, 0, c, NULL,
                                interpolate_partition_fn, iquery_points);

#ifdef FCLAW_ENABLE_MPI
    /* notify the producer about the point-array-messages it will receive,
     * allocate an receive buffer according to the transmitted payloads and
     * post Irecvs for the point-arrays */
    consumer_producer_notify (g);

    /* post Isends for the point-arrays as well as Irecvs for the updated
     * point-arrays containing the interpolation prodata */
    consumer_post_messages (c);

    fclaw_global_essentialf ("OVERLAP: producer local search\n");

    /* interpolate the point-arrays as soon as they arrive and send them back to
     * the consumer side in a non-blocking way */
    producer_interpolate (p);

    fclaw_global_essentialf ("OVERLAP: consumer query point update\n");

    /* compute the interpolation data of the query points based on the
     * updated point-arrays */
    consumer_update_query_points (c);

    /* wait for the communication to complete */
    consumer_waitall (c);
    producer_waitall (p);
#else
    c->iconrank = 0;
#endif

    fclaw_global_essentialf ("OVERLAP: local interpolation\n");

    /* local, in-place part of the interpolation */
    consumer_producer_update_local (g);

    /* free remaining communication data */
    consumer_free_communication_data (c);
    producer_free_communication_data (p);
    sc_array_destroy (iquery_points);
}
