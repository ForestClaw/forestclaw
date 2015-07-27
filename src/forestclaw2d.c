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

#include <forestclaw2d.h>
#include <p4est_bits.h>
#include <p4est_wrap.h>

#define FCLAW2D_DOMAIN_TAG_SERIALIZE 4526

const fclaw2d_patch_flags_t fclaw2d_patch_block_face_flags[4] = {
    FCLAW2D_PATCH_ON_BLOCK_FACE_0,
    FCLAW2D_PATCH_ON_BLOCK_FACE_1,
    FCLAW2D_PATCH_ON_BLOCK_FACE_2,
    FCLAW2D_PATCH_ON_BLOCK_FACE_3
};

/* This is already deprecated:
 * I need to go back to fclaw_base and see how to make that usable. */
void
fclaw2d_global_log (int log_priority, const char *message)
{
    /* TODO: establish an fclaw_package_id */
    SC_GEN_LOG (sc_package_id, SC_LC_GLOBAL, log_priority, message);
}

double
fclaw2d_domain_global_maximum (fclaw2d_domain_t * domain, double d)
{
    int mpiret;
    double gd;

    mpiret = sc_MPI_Allreduce (&d, &gd, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
                               domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return gd;
}

double
fclaw2d_domain_global_sum (fclaw2d_domain_t * domain, double d)
{
    int mpiret;
    double gd;

    mpiret = sc_MPI_Allreduce (&d, &gd, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                               domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return gd;
}

void
fclaw2d_domain_barrier (fclaw2d_domain_t * domain)
{
    int mpiret;

    mpiret = sc_MPI_Barrier (domain->mpicomm);
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
    FCLAW_ASSERT (0 <= icorner && icorner < P4EST_CHILDREN);
    faces[0] = p4est_corner_faces[icorner][0];
    faces[1] = p4est_corner_faces[icorner][1];
}

int
fclaw2d_patch_corner_dimension (const fclaw2d_patch_t * patch, int cornerno)
{
    const int childid = fclaw2d_patch_childid (patch);

    FCLAW_ASSERT (0 <= cornerno && cornerno < P4EST_CHILDREN);

    return (patch->level == 0 ||
            cornerno == childid ||
            cornerno == P4EST_CHILDREN - 1 - childid) ? 0 : 1;
}

int
fclaw2d_patch_childid (const fclaw2d_patch_t * patch)
{
    const int childid = patch->flags & FCLAW2D_PATCH_CHILDID;

    FCLAW_ASSERT (0 <= childid && childid < P4EST_CHILDREN);

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
fclaw2d_domain_attribute_add (fclaw2d_domain_t * domain,
                              const char *name, void *attribute)
{
    sc_keyvalue_t *a = domain->attributes;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (!sc_keyvalue_exists (a, name));
    sc_keyvalue_set_pointer (a, name, attribute);
}

void *
fclaw2d_domain_attribute_access (fclaw2d_domain_t * domain,
                                 const char *name, void *default_attr)
{
    sc_keyvalue_t *a = domain->attributes;

    FCLAW_ASSERT (a != NULL);
    return sc_keyvalue_get_pointer (a, name, default_attr);
}

void
fclaw2d_domain_attribute_remove (fclaw2d_domain_t * domain, const char *name)
{
    sc_keyvalue_t *a = domain->attributes;
#ifdef FCLAW_ENABLE_DEBUG
    sc_keyvalue_entry_type_t et;
#endif

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (sc_keyvalue_exists (a, name));
#ifndef FCLAW_ENABLE_DEBUG
    (void)
#else
    et =
#endif
        sc_keyvalue_unset (a, name);
    FCLAW_ASSERT (et == SC_KEYVALUE_ENTRY_POINTER);
}

fclaw2d_patch_t *
fclaw2d_domain_get_patch (fclaw2d_domain_t * domain, int blockno, int patchno)
{
    fclaw2d_block_t *block;

    /* remote patch */
    if (blockno == -1)
    {
        FCLAW_ASSERT (0 <= patchno && patchno < domain->num_ghost_patches);
        return domain->ghost_patches + patchno;
    }

    /* local patch */
    FCLAW_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    block = domain->blocks + blockno;

    FCLAW_ASSERT (0 <= patchno && patchno < block->num_patches);
    return block->patches + patchno;
}

void
fclaw2d_domain_iterate_level (fclaw2d_domain_t * domain, int level,
                              fclaw2d_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    FCLAW_ASSERT (0 <= level && level <= domain->possible_maxlevel);
    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        for (patch = block->patchbylevel[level];
             patch != NULL; patch = patch->u.next)
        {
            j = (int) (patch - block->patches);
            FCLAW_ASSERT (0 <= j && j < block->num_patches);
            FCLAW_ASSERT (patch->level == level);
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
#ifdef FCLAW_ENABLE_DEBUG
                int k;
                for (k = 0; k < P4EST_CHILDREN; ++k)
                {
                    FCLAW_ASSERT (j + k < block->num_patches);
                    FCLAW_ASSERT (fclaw2d_patch_childid (patch + k) == k);
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
#ifdef FCLAW_ENABLE_DEBUG
    fclaw2d_block_t *block;
#endif

    FCLAW_ASSERT (domain->pp_owned);
    anyboundary = 0;

    FCLAW_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    FCLAW_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    FCLAW_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

#ifdef FCLAW_ENABLE_DEBUG
    block = domain->blocks + blockno;
#endif
    FCLAW_ASSERT (0 <= patchno && patchno < block->num_patches);

    /* TODO: Use patch block face flags as shortcut */

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    totalleaf = tree->quadrants_offset + (p4est_locidx_t) patchno;
    FCLAW_ASSERT (0 <= totalleaf && totalleaf < p4est->local_num_quadrants);
    for (faceno = 0; faceno < P4EST_FACES; ++faceno)
    {
        qtq = mesh->quad_to_quad[P4EST_FACES * totalleaf + faceno];
        qtf = mesh->quad_to_face[P4EST_FACES * totalleaf + faceno];
        if (qtq == totalleaf && qtf == faceno)
        {
            FCLAW_ASSERT (block->is_boundary[faceno]);
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

    FCLAW_ASSERT (domain->pp_owned);
    FCLAW_ASSERT (0 <= qtq);
    FCLAW_ASSERT (qtq <
                  mesh->local_num_quadrants + mesh->ghost_num_quadrants);
    if (qtq < mesh->local_num_quadrants)
    {
        /* processor-local neighbor */
        *proc = domain->mpirank;
        *blockno = (int) mesh->quad_to_tree[qtq];
        FCLAW_ASSERT ((int) wrap->p4est->first_local_tree <= *blockno);
        FCLAW_ASSERT (*blockno <= (int) wrap->p4est->last_local_tree);
        block = domain->blocks + *blockno;
        qtq -= block->num_patches_before;
        FCLAW_ASSERT (0 <= qtq && qtq < block->num_patches);
        *patchno = (int) qtq;   /* patch number within the block as usual */
    }
    else
    {
        /* off-processor ghost neighbor */
        qtq -= mesh->local_num_quadrants;
        FCLAW_ASSERT (qtq >= 0 && qtq < mesh->ghost_num_quadrants);
        *proc = mesh->ghost_to_proc[qtq];
        gq = p4est_quadrant_array_index (&ghost->ghosts, qtq);
        FCLAW_ASSERT (0 <= gq->p.piggy3.which_tree);
        FCLAW_ASSERT (gq->p.piggy3.which_tree < wrap->conn->num_trees);
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
#ifdef FCLAW_ENABLE_DEBUG
    fclaw2d_block_t *block;
#endif

    FCLAW_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);

    FCLAW_ASSERT (domain->pp_owned);

    FCLAW_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    FCLAW_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    FCLAW_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

#ifdef FCLAW_ENABLE_DEBUG
    block = domain->blocks + blockno;
#endif
    FCLAW_ASSERT (0 <= patchno && patchno < block->num_patches);

    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    totalleaf = tree->quadrants_offset + (p4est_locidx_t) patchno;
    FCLAW_ASSERT (0 <= totalleaf && totalleaf < p4est->local_num_quadrants);
    FCLAW_ASSERT (0 <= faceno && faceno < P4EST_FACES);

    qtq = mesh->quad_to_quad[P4EST_FACES * totalleaf + faceno];
    qtf = mesh->quad_to_face[P4EST_FACES * totalleaf + faceno];
    if (qtq == totalleaf && qtf == faceno)
    {
        /* physical domain boundary encoded by same patch face */
        FCLAW_ASSERT (block->is_boundary[faceno]);
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
            FCLAW_ASSERT (k == 0 || hblockno[k - 1] == hblockno[k]);
        }
        *rblockno = hblockno[0];
        *rfaceno = qtf + num_orientations;
        FCLAW_ASSERT (*rfaceno >= 0);
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
            FCLAW_ASSERT (0 <= *rfaceno && *rfaceno < num_orientations);
            return FCLAW2D_PATCH_SAMESIZE;
        }
        else
        {
            /* double-size neighbor */
            *rfaceno = (int) qtf % num_orientations;
            /* the number of our patch within the bigger neighbor subfaces */
            rproc[1] = (int) qtf / num_orientations - 1;
            FCLAW_ASSERT (0 <= rproc[1] && rproc[1] < P4EST_HALF);
            return FCLAW2D_PATCH_DOUBLESIZE;
        }
    }
}

void
fclaw2d_patch_face_swap (int *faceno, int *rfaceno)
{
    const int iface = *faceno;
    const int orientation = *rfaceno / P4EST_FACES;

    P4EST_ASSERT (orientation < P4EST_HALF);
    P4EST_ASSERT (*faceno < P4EST_FACES);

    *faceno = *rfaceno % P4EST_FACES;
    *rfaceno = iface + P4EST_FACES * orientation;

    P4EST_ASSERT (*faceno < P4EST_FACES);
}

void
fclaw2d_patch_face_transformation (int faceno, int rfaceno, int ftransform[])
{
    p4est_expand_face_transform (faceno, rfaceno, ftransform);
}

void
fclaw2d_patch_transform_face (fclaw2d_patch_t * ipatch,
                              fclaw2d_patch_t * opatch,
                              const int ftransform[],
                              int mx, int my, int based, int *i, int *j)
{
    double Rmx;

    FCLAW_ASSERT (ipatch->level == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && ipatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);

    FCLAW_ASSERT (mx >= 1 && mx == my);
    FCLAW_ASSERT (based == 0 || based == 1);

#if 0
    printf ("Test I: IP %g %g %d FT %d %d %d %d %d %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level,
            ftransform[0], ftransform[2], ftransform[3], ftransform[5],
            ftransform[6], ftransform[8], mx, *i, *j, based);
#endif

    /* work with doubles -- exact for integers up to 52 bits of precision */
    Rmx = (double) mx *(double) (1 << ipatch->level);

    if (ftransform[8] & 4)
    {
        /* The two patches are in the same block.  ftransform is not used */
        *i += (int) ((ipatch->xlower - opatch->xlower) * Rmx);
        *j += (int) ((ipatch->ylower - opatch->ylower) * Rmx);
    }
    else
    {
        const int *my_axis = &ftransform[0];
        const int *target_axis = &ftransform[3];
        const int *edge_reverse = &ftransform[6];
        double my_xyz[2], target_xyz[2];

        /* the reference cube is stretched to mx times my units */
        my_xyz[0] = ipatch->xlower * Rmx + *i - .5 * based;
        my_xyz[1] = ipatch->ylower * Rmx + *j - .5 * based;

        /* transform transversal direction */
        target_xyz[target_axis[0]] =
            !edge_reverse[0] ? my_xyz[my_axis[0]] : Rmx - my_xyz[my_axis[0]];

        /* transform normal direction */
        switch (edge_reverse[2])
        {
        case 0:
            target_xyz[target_axis[2]] = -my_xyz[my_axis[2]];
            break;
        case 1:
            target_xyz[target_axis[2]] = my_xyz[my_axis[2]] + Rmx;
            break;
        case 2:
            target_xyz[target_axis[2]] = my_xyz[my_axis[2]] - Rmx;
            break;
        case 3:
            target_xyz[target_axis[2]] = 2. * Rmx - my_xyz[my_axis[2]];
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* transform back to integer coordinates: this is exact */
        *i = (int) (target_xyz[0] - opatch->xlower * Rmx + .5 * based);
        *j = (int) (target_xyz[1] - opatch->ylower * Rmx + .5 * based);
    }

#if 0
    printf ("Test O: IP %g %g IJ %d %d\n",
            opatch->xlower, opatch->ylower, *i, *j);
#endif
}

void
fclaw2d_patch_transform_face2 (fclaw2d_patch_t * ipatch,
                               fclaw2d_patch_t * opatch,
                               const int ftransform[],
                               int mx, int my, int based, int i[], int j[])
{
    int kt, kn;
    int di, dj;
    double Rmx;

    FCLAW_ASSERT (ipatch->level + 1 == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && opatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);

    FCLAW_ASSERT (mx >= 1 && mx == my);
    FCLAW_ASSERT (based == 0 || based == 1);

#if 0
    printf ("Test I: IP %g %g %d FT %d %d %d %d %d %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level,
            ftransform[0], ftransform[2], ftransform[3], ftransform[5],
            ftransform[6], ftransform[8], mx, *i, *j, based);
#endif

    /* work with doubles -- exact for integers up to 52 bits of precision */
    Rmx = (double) mx *(double) (1 << opatch->level);

    if (ftransform[8] & 4)
    {
        /* The two patches are in the same block.  ftransform is undefined */
        di = based + (int)
            ((ipatch->xlower - opatch->xlower) * Rmx + 2. * (*i - based));
        dj = based + (int)
            ((ipatch->ylower - opatch->ylower) * Rmx + 2. * (*j - based));

        /* In the same block, the order of child cells is canonical */
        for (kt = 0; kt < 2; ++kt)
        {
            for (kn = 0; kn < 2; ++kn)
            {
                i[2 * kt + kn] = di + kn;
                j[2 * kt + kn] = dj + kt;
            }
        }
    }
    else
    {
        const int *my_axis = &ftransform[0];
        const int *target_axis = &ftransform[3];
        const int *edge_reverse = &ftransform[6];
        int bt, bn;
        int nx, ny;
        int sx, sy;
        double my_xyz[2], target_xyz[2];

        /* the reference cube is stretched to mx times my units */
        my_xyz[0] = ipatch->xlower * Rmx + 2. * (*i + .5 - based);
        my_xyz[1] = ipatch->ylower * Rmx + 2. * (*j + .5 - based);

        /* transform transversal direction */
        target_xyz[target_axis[0]] =
            !edge_reverse[0] ? my_xyz[my_axis[0]] : Rmx - my_xyz[my_axis[0]];

        /* transform normal direction */
        switch (edge_reverse[2])
        {
        case 0:
            target_xyz[target_axis[2]] = -my_xyz[my_axis[2]];
            break;
        case 1:
            target_xyz[target_axis[2]] = my_xyz[my_axis[2]] + Rmx;
            break;
        case 2:
            target_xyz[target_axis[2]] = my_xyz[my_axis[2]] - Rmx;
            break;
        case 3:
            target_xyz[target_axis[2]] = 2. * Rmx - my_xyz[my_axis[2]];
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* move back into integer coordinates: this is exact */
        di = (int) (target_xyz[0] - opatch->xlower * Rmx) + based - 1;
        dj = (int) (target_xyz[1] - opatch->ylower * Rmx) + based - 1;

        /* Run through the child cells in order of the small patch */
        for (sy = 0; sy < 2; ++sy)
        {
            for (sx = 0; sx < 2; ++sx)
            {
                /* Compute small patch coordinate in (transverse, normal) order */
                if (target_axis[0])
                {
                    kt = sy;
                    kn = sx;
                }
                else
                {
                    kt = sx;
                    kn = sy;
                }

                /* Transform into big patch (transverse, normal) coordinate */
                bt = !edge_reverse[0] ? kt : !kt;
                bn = (edge_reverse[2] == 1
                      || edge_reverse[2] == 2) ? kn : !kn;

                /* Compute coordinate relative to the big patch */
                if (my_axis[0])
                {
                    nx = bn;
                    ny = bt;
                }
                else
                {
                    nx = bt;
                    ny = bn;
                }

                /* assign value in proper place */
                i[2 * ny + nx] = di + sx;
                j[2 * ny + nx] = dj + sy;
            }
        }
    }

#if 0
    printf ("Test O: OP %g %g %d I %d %d %d %d J %d %d %d %d\n",
            opatch->xlower, opatch->ylower, opatch->level,
            i[0], i[1], i[2], i[3], j[0], j[1], j[2], j[3]);
#endif
}

static inline p4est_locidx_t
fclaw2d_array_index_locidx (sc_array_t * array, p4est_locidx_t li)
{
    return *(p4est_locidx_t *) sc_array_index_int (array, (int) li);
}

int
fclaw2d_patch_corner_neighbors (fclaw2d_domain_t * domain,
                                int blockno, int patchno, int cornerno,
                                int *rproc, int *rblockno, int *rpatchno,
                                int *rcorner,
                                fclaw2d_patch_relation_t * neighbor_size)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    p4est_locidx_t local_num, qid;
    p4est_locidx_t cornerid, cstart, cend;
    const p4est_quadrant_t *q;
    p4est_tree_t *rtree;
    fclaw2d_block_t *block;

    FCLAW_ASSERT (domain->num_ghost_patches ==
                  (int) mesh->ghost_num_quadrants);

    FCLAW_ASSERT (domain->pp_owned);

    FCLAW_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    FCLAW_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    FCLAW_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);

    block = domain->blocks + blockno;
    FCLAW_ASSERT (0 <= patchno && patchno < block->num_patches);
    FCLAW_ASSERT (0 <= cornerno && cornerno < P4EST_CHILDREN);

    local_num = block->num_patches_before + patchno;
    qid = mesh->quad_to_corner[P4EST_CHILDREN * local_num + cornerno];

    /* We are not yet ready for general multiblock connectivities where more
     * than four blocks meet at a corner */
    if (qid >= 0)
    {
        FCLAW_ASSERT (0 <= qid);
        if (qid >= mesh->local_num_quadrants + mesh->ghost_num_quadrants)
        {
            /* This is an inter-tree (face or corner) corner neighbor */
            cornerid =
                qid - (mesh->local_num_quadrants + mesh->ghost_num_quadrants);
            FCLAW_ASSERT (cornerid < mesh->local_num_corners);
            cstart =
                fclaw2d_array_index_locidx (mesh->corner_offset, cornerid);
            cend =
                fclaw2d_array_index_locidx (mesh->corner_offset,
                                            cornerid + 1);
            if (cstart + 1 < cend)
            {
                /* At least a five-corner, which is currently not supported */
                qid = -1;
            }
            else
            {
                FCLAW_ASSERT (cstart + 1 == cend);
                qid = fclaw2d_array_index_locidx (mesh->corner_quad, cstart);
                *rcorner = (int)
                    *(int8_t *) sc_array_index_int (mesh->corner_corner,
                                                    (int) cstart);
                FCLAW_ASSERT (0 <= *rcorner && *rcorner < P4EST_CHILDREN);
            }
        }
        else
        {
            /* for intra-tree corners we take the corner is opposite */
            *rcorner = cornerno ^ (P4EST_CHILDREN - 1);
        }
    }

    if (qid < 0)
    {
        /* The value -1 is expected for a corner on the physical boundary */
        /* Currently we also return this for five- and more-corners */
        *neighbor_size = FCLAW2D_PATCH_BOUNDARY;
        *rcorner = -1;
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
            *neighbor_size = FCLAW2D_PATCH_DOUBLESIZE;
            break;
        case 0:
            *neighbor_size = FCLAW2D_PATCH_SAMESIZE;
            break;
        case 1:
            *neighbor_size = FCLAW2D_PATCH_HALFSIZE;
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

    return *neighbor_size != FCLAW2D_PATCH_BOUNDARY;
}

void
fclaw2d_patch_corner_swap (int *cornerno, int *rcornerno)
{
    int swap;

    swap = *cornerno;
    *cornerno = *rcornerno;
    *rcornerno = swap;
}

void
fclaw2d_patch_transform_corner (fclaw2d_patch_t * ipatch,
                                fclaw2d_patch_t * opatch,
                                int icorner, int is_block_boundary,
                                int mx, int my, int based, int *i, int *j)
{
    double Rmx, xshift, yshift;

    FCLAW_ASSERT (ipatch->level == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && ipatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);

    FCLAW_ASSERT (mx >= 1 && mx == my);
    FCLAW_ASSERT (based == 0 || based == 1);

#if 0
    printf ("Test I: IP %g %g %d IC %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level, icorner,
            mx, *i, *j, based);
#endif

    /* Work with doubles -- exact for integers up to 52 bits of precision */
    Rmx = (double) mx *(double) (1 << ipatch->level);
    if (!is_block_boundary)
    {
        /* The lower left coordinates are with respect to the same origin */
        xshift = yshift = 0.;
    }
    else
    {
        /* We need to add/substract the shift due to translation of one block */
        if ((icorner & 1) == 0)
        {
            /* This corner is on the left face of the patch */
            xshift = +1.;
        }
        else
        {
            /* This corner is on the right face of the patch */
            xshift = -1.;
        }
        if ((icorner & 2) == 0)
        {
            /* This corner is on the bottom face of the patch */
            yshift = +1.;
        }
        else
        {
            /* This corner is on the top face of the patch */
            yshift = -1.;
        }
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    *i += (int) ((ipatch->xlower - opatch->xlower + xshift) * Rmx);
    *j += (int) ((ipatch->ylower - opatch->ylower + yshift) * Rmx);
}

void
fclaw2d_patch_transform_corner2 (fclaw2d_patch_t * ipatch,
                                 fclaw2d_patch_t * opatch,
                                 int icorner, int is_block_boundary,
                                 int mx, int my, int based, int i[], int j[])
{
    int kt, kn;
    int di, dj;
    double Rmx, xshift, yshift;

    FCLAW_ASSERT (ipatch->level + 1 == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && opatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);

    FCLAW_ASSERT (mx >= 1 && mx == my);
    FCLAW_ASSERT (based == 0 || based == 1);

#if 0
    printf ("Test I: IP %g %g %d IC %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level, icorner,
            mx, *i, *j, based);
#endif

    /* work with doubles -- exact for integers up to 52 bits of precision */
    Rmx = (double) mx *(double) (1 << opatch->level);
    if (!is_block_boundary)
    {
        /* The lower left coordinates are with respect to the same origin */
        xshift = yshift = 0.;
    }
    else
    {
        /* We need to add/substract the shift due to translation of one block */
        if ((icorner & 1) == 0)
        {
            /* This corner is on the left face of the patch */
            xshift = +1.;
        }
        else
        {
            /* This corner is on the right face of the patch */
            xshift = -1.;
        }
        if ((icorner & 2) == 0)
        {
            /* This corner is on the bottom face of the patch */
            yshift = +1.;
        }
        else
        {
            /* This corner is on the top face of the patch */
            yshift = -1.;
        }
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    di = based
        + (int) ((ipatch->xlower - opatch->xlower + xshift) * Rmx +
                 2. * (*i - based));
    dj = based
        + (int) ((ipatch->ylower - opatch->ylower + yshift) * Rmx +
                 2. * (*j - based));

    /* Without any rotation, the order of child cells is canonical */
    for (kt = 0; kt < 2; ++kt)
    {
        for (kn = 0; kn < 2; ++kn)
        {
            i[2 * kt + kn] = di + kn;
            j[2 * kt + kn] = dj + kt;
        }
    }
}

void
fclaw2d_domain_set_refinement (fclaw2d_domain_t * domain,
                               int smooth_refine, int coarsen_delay)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (coarsen_delay >= 0);

    domain->smooth_refine = smooth_refine;
    p4est_wrap_set_coarsen_delay (wrap, coarsen_delay, 0);
}

void
fclaw2d_patch_mark_refine (fclaw2d_domain_t * domain, int blockno,
                           int patchno)
{
    fclaw2d_patch_t *patch;
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    /* bump target level up */
    patch = fclaw2d_domain_get_patch (domain, blockno, patchno);
    patch->target_level = patch->level + 1;
    patch->target_level = SC_MIN (patch->target_level, P4EST_QMAXLEVEL);

    /* if we do smooth refinement, all marking is done inside adapt */
    if (!domain->smooth_refine)
    {
        p4est_wrap_mark_refine (wrap,
                                (p4est_locidx_t) blockno,
                                (p4est_locidx_t) patchno);
    }
}

void
fclaw2d_patch_mark_coarsen (fclaw2d_domain_t * domain, int blockno,
                            int patchno)
{
    fclaw2d_patch_t *patch;
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    /* bump target level down */
    patch = fclaw2d_domain_get_patch (domain, blockno, patchno);
    patch->target_level = patch->level - 1;
    patch->target_level = SC_MAX (patch->target_level, 0);

    /* if we do smooth refinement, all marking is done inside adapt */
    if (!domain->smooth_refine)
    {
        p4est_wrap_mark_coarsen (wrap,
                                 (p4est_locidx_t) blockno,
                                 (p4est_locidx_t) patchno);
    }
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

    FCLAW_ASSERT (!old_domain->pp_owned);
    FCLAW_ASSERT (new_domain->pp_owned);
    FCLAW_ASSERT (old_domain->pp == new_domain->pp);
    FCLAW_ASSERT (old_domain->num_blocks == new_domain->num_blocks);
    for (i = 0; i < old_domain->num_blocks; i++)
    {
        old_block = old_domain->blocks + i;
        new_block = new_domain->blocks + i;
        for (oj = nj = 0; oj < old_block->num_patches;)
        {
            FCLAW_ASSERT (nj < new_block->num_patches);
            old_patch = old_block->patches + oj;
            new_patch = new_block->patches + nj;
            FCLAW_ASSERT (abs (old_patch->level - new_patch->level) <= 1);
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
        FCLAW_ASSERT (oj == old_block->num_patches);
        FCLAW_ASSERT (nj == new_block->num_patches);
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
            FCLAW_ASSERT (zz ==
                          (size_t) (block->num_patches_before + patchno));

            q = p4est_quadrant_array_index (&tree->quadrants,
                                            (p4est_locidx_t) patchno);
            patch_data[zz] = q->p.user_data;
        }
    }
    FCLAW_ASSERT (zz == (size_t) domain->local_num_patches);
}

void
fclaw2d_domain_allocate_before_partition (fclaw2d_domain_t * domain,
                                          size_t data_size,
                                          void ***patch_data)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (*patch_data == NULL);

    p4est_reset_data (wrap->p4est, data_size, NULL,
                      wrap->p4est->user_pointer);

    *patch_data = FCLAW_ALLOC (void *, domain->local_num_patches);
    fclaw2d_domain_assign_for_partition (domain, *patch_data);
}

void
fclaw2d_domain_retrieve_after_partition (fclaw2d_domain_t * domain,
                                         void ***patch_data)
{
    *patch_data =
        FCLAW_REALLOC (*patch_data, void *, domain->local_num_patches);
    fclaw2d_domain_assign_for_partition (domain, *patch_data);
}

void
fclaw2d_domain_iterate_partitioned (fclaw2d_domain_t * old_domain,
                                    fclaw2d_domain_t * new_domain,
                                    fclaw2d_transfer_callback_t tcb,
                                    void *user)
{
    int uf, ul, uof;
    int obp, nbp, oskip, nskip;
    int opre, npre, ostay, nstay;
    int i, oj, nj;
#ifdef FCLAW_ENABLE_DEBUG
    int odone, ndone;
#endif
    fclaw2d_block_t *old_block, *new_block;
    fclaw2d_patch_t *old_patch, *new_patch;

    FCLAW_ASSERT (!old_domain->pp_owned);
    FCLAW_ASSERT (new_domain->pp_owned);
    FCLAW_ASSERT (new_domain->just_partitioned);
    FCLAW_ASSERT (old_domain->pp == new_domain->pp);
    FCLAW_ASSERT (old_domain->num_blocks == new_domain->num_blocks);

    /* get information on overlapping window of local patches */
    uf = new_domain->partition_unchanged_first;
    ul = new_domain->partition_unchanged_length;
    uof = new_domain->partition_unchanged_old_first;

    /* go through the blocks in the old and new partitions */
    for (i = 0; i < new_domain->num_blocks; i++)
    {
        /* if this block has no patches we skip it */
        new_block = new_domain->blocks + i;
        nbp = new_block->num_patches;
        if (nbp == 0)
        {
            continue;
        }
        nskip = new_block->num_patches_before;
#ifdef FCLAW_ENABLE_DEBUG
        ndone = 0;
#endif

        /* This block has patches.  See which of these preexisted locally */
        npre = uf - nskip;
        npre = SC_MAX (npre, 0);
        if (npre >= nbp)
        {
            npre = nstay = nbp;
        }
        else
        {
            nstay = uf + ul - nskip;
            nstay = SC_MAX (nstay, 0);
            nstay = SC_MIN (nstay, nbp);
        }
        FCLAW_ASSERT (0 <= npre && npre <= nstay && nstay <= nbp);

        /* Iterate over first set of newly received patches */
        for (nj = 0; nj < npre; ++nj)
        {
            FCLAW_ASSERT (nj < new_block->num_patches);
            new_patch = new_block->patches + nj;
            tcb (old_domain, NULL, new_domain, new_patch, i, -1, nj, user);
#ifdef FCLAW_ENABLE_DEBUG
            ++ndone;
#endif
        }

        /* are we done with the patches in this block, all of which new? */
        if (nj == nbp)
        {
            continue;
        }
        old_block = old_domain->blocks + i;
        obp = old_block->num_patches;
        oskip = old_block->num_patches_before;
#ifdef FCLAW_ENABLE_DEBUG
        odone = 0;
#endif

        /* We may have preexisting patches.  Figure out the source range */
        opre = uof - oskip;
        opre = SC_MAX (opre, 0);
        if (opre >= obp)
        {
            opre = ostay = obp;
        }
        else
        {
            ostay = uof + ul - oskip;
            ostay = SC_MAX (ostay, 0);
            ostay = SC_MIN (ostay, obp);
        }
        FCLAW_ASSERT (0 <= opre && opre <= ostay && ostay <= obp);

        /* Iterate over the patches that stay local, if any */
        FCLAW_ASSERT (ostay - opre == nstay - npre);
        for (oj = opre; nj < nstay; ++oj, ++nj)
        {
            FCLAW_ASSERT (oj < old_block->num_patches);
            FCLAW_ASSERT (nj < new_block->num_patches);
            old_patch = old_block->patches + oj;
            new_patch = new_block->patches + nj;
            tcb (old_domain, old_patch, new_domain, new_patch,
                 i, oj, nj, user);
#ifdef FCLAW_ENABLE_DEBUG
            ++odone;
            ++ndone;
#endif
        }
        FCLAW_ASSERT (oj == ostay);
        FCLAW_ASSERT (odone <= obp);

        /* Iterate over second set of newly received patches */
        for (; nj < nbp; ++nj)
        {
            FCLAW_ASSERT (nj < new_block->num_patches);
            new_patch = new_block->patches + nj;
            tcb (old_domain, NULL, new_domain, new_patch, i, -1, nj, user);
#ifdef FCLAW_ENABLE_DEBUG
            ++ndone;
#endif
        }

        /* consistency checks */
        FCLAW_ASSERT (nj == new_block->num_patches);
        FCLAW_ASSERT (ndone == nbp);
    }
}

void
fclaw2d_domain_free_after_partition (fclaw2d_domain_t * domain,
                                     void ***patch_data)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_FREE (*patch_data);
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

    e = FCLAW_ALLOC (fclaw2d_domain_exchange_t, 1);
    e->data_size = data_size;
    e->num_exchange_patches = domain->num_exchange_patches;
    e->num_ghost_patches = domain->num_ghost_patches;

    e->patch_data = FCLAW_ALLOC (void *, domain->num_exchange_patches);
    e->ghost_data = FCLAW_ALLOC (void *, domain->num_ghost_patches);
    e->ghost_contiguous_memory = m = FCLAW_ALLOC (char,
                                                  (size_t)
                                                  domain->num_ghost_patches *
                                                  data_size);
    for (i = 0; i < domain->num_ghost_patches; ++i)
    {
        e->ghost_data[i] = m;
        m += data_size;
    }

    e->async_state = NULL;
    e->by_levels = 0;
    e->inside_async = 0;

    return e;
}

void
fclaw2d_domain_ghost_exchange (fclaw2d_domain_t * domain,
                               fclaw2d_domain_exchange_t * e,
                               int exchange_minlevel, int exchange_maxlevel)
{
    FCLAW_ASSERT (e->async_state == NULL);
    FCLAW_ASSERT (!e->inside_async);

    fclaw2d_domain_ghost_exchange_begin
        (domain, e, exchange_minlevel, exchange_maxlevel);
    fclaw2d_domain_ghost_exchange_end (domain, e);
}

void
fclaw2d_domain_ghost_exchange_begin (fclaw2d_domain_t * domain,
                                     fclaw2d_domain_exchange_t * e,
                                     int exchange_minlevel,
                                     int exchange_maxlevel)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_ghost_t *ghost = wrap->match_aux ? wrap->ghost_aux : wrap->ghost;
    p4est_ghost_exchange_t *exc;

    /* we must not be in an active exchange already */
    FCLAW_ASSERT (e->async_state == NULL);
    FCLAW_ASSERT (!e->inside_async);

#if 0
    P4EST_LDEBUGF ("Patches exchange begin %d %d, ghost %d %d\n",
                   e->num_exchange_patches, (int) ghost->mirrors.elem_count,
                   e->num_ghost_patches, (int) ghost->ghosts.elem_count);
#endif

    FCLAW_ASSERT (e->num_exchange_patches == (int) ghost->mirrors.elem_count);
    FCLAW_ASSERT (e->num_ghost_patches == (int) ghost->ghosts.elem_count);
    if (exchange_minlevel <= domain->global_minlevel &&
        domain->global_maxlevel <= exchange_maxlevel)
    {
        e->by_levels = 0;
        exc = p4est_ghost_exchange_custom_begin
            (wrap->p4est, ghost,
             e->data_size, e->patch_data, e->ghost_contiguous_memory);
    }
    else
    {
        e->by_levels = 1;
        exc = p4est_ghost_exchange_custom_levels_begin
            (wrap->p4est, ghost,
             exchange_minlevel, exchange_maxlevel, e->data_size,
             e->patch_data, e->ghost_contiguous_memory);
    }
    e->async_state = exc;
    e->inside_async = 1;
}

void
fclaw2d_domain_ghost_exchange_end (fclaw2d_domain_t * domain,
                                   fclaw2d_domain_exchange_t * e)
{
    p4est_ghost_exchange_t *exc = (p4est_ghost_exchange_t *) e->async_state;

    FCLAW_ASSERT (exc != NULL);
    FCLAW_ASSERT (e->inside_async);

    if (!e->by_levels)
    {
        p4est_ghost_exchange_custom_end (exc);
    }
    else
    {
        p4est_ghost_exchange_custom_levels_end (exc);
    }

    e->async_state = NULL;
    e->inside_async = 0;
}

void
fclaw2d_domain_free_after_exchange (fclaw2d_domain_t * domain,
                                    fclaw2d_domain_exchange_t * e)
{
    FCLAW_FREE (e->ghost_contiguous_memory);
    FCLAW_FREE (e->ghost_data);
    FCLAW_FREE (e->patch_data);
    FCLAW_FREE (e);
}

void
fclaw2d_domain_serialization_enter (fclaw2d_domain_t * domain)
{
    int mpiret;
    int i;
    sc_MPI_Status status;

    if (domain->mpirank > 0)
    {
        mpiret = sc_MPI_Recv (&i, 1, sc_MPI_INT, domain->mpirank - 1,
                              FCLAW2D_DOMAIN_TAG_SERIALIZE, domain->mpicomm,
                              &status);
        SC_CHECK_MPI (mpiret);
        FCLAW_ASSERT (i == 0);
    }
}

void
fclaw2d_domain_serialization_leave (fclaw2d_domain_t * domain)
{
    int mpiret;
    int i = 0;

    if (domain->mpirank + 1 < domain->mpisize)
    {
        mpiret = sc_MPI_Send (&i, 1, sc_MPI_INT, domain->mpirank + 1,
                              FCLAW2D_DOMAIN_TAG_SERIALIZE, domain->mpicomm);
        SC_CHECK_MPI (mpiret);
    }
}
