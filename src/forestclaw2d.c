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

#ifndef P4_TO_P8
#include <forestclaw2d.h>
#include <p4est_bits.h>
#include <p4est_wrap.h>
#else
#include <forestclaw3d.h>
#include <p8est_bits.h>
#include <p8est_wrap.h>
#endif

#ifndef P4_TO_P8
#define FCLAW2D_DOMAIN_TAG_SERIALIZE 4526

const fclaw2d_patch_flags_t fclaw2d_patch_block_face_flags[4] = {
    FCLAW2D_PATCH_ON_BLOCK_FACE_0,
    FCLAW2D_PATCH_ON_BLOCK_FACE_1,
    FCLAW2D_PATCH_ON_BLOCK_FACE_2,
    FCLAW2D_PATCH_ON_BLOCK_FACE_3
};

#else
#define FCLAW2D_DOMAIN_TAG_SERIALIZE 4527
#endif

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
                             int icorner, int faces[P4EST_DIM])
{
    FCLAW_ASSERT (0 <= icorner && icorner < P4EST_CHILDREN);
    faces[0] = p4est_corner_faces[icorner][0];
    faces[1] = p4est_corner_faces[icorner][1];
#ifdef P4_TO_P8
    faces[2] = p8est_corner_faces[icorner][2];
#endif
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

static fclaw2d_patch_t *
fclaw2d_domain_get_patch (fclaw2d_domain_t * domain, int blockno, int patchno)
{
    fclaw2d_block_t *block;

#if 0
    /* remote patch */
    if (blockno == -1)
    {
        FCLAW_ASSERT (0 <= patchno && patchno < domain->num_ghost_patches);
        return domain->ghost_patches + patchno;
    }
#endif

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

#ifdef FCLAW_ENABLE_DEBUG
#ifndef P4_TO_P8
static const int normal_out[4] = { 0, 1, 0, 1 };
#else
static const int normal_out[6] = { 0, 1, 0, 1, 0, 1 };
#endif
#endif

int
fclaw2d_patch_normal_match (fclaw2d_domain_t * domain,
                            int blockno, int patchno, int faceno)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;
    int qtfi;
    p4est_locidx_t totalleaf;
    p4est_tree_t *tree;
#ifdef FCLAW_ENABLE_DEBUG
    fclaw2d_block_t *block;
    int num_orient = fclaw2d_domain_num_orientations (domain);
#endif

    /* are we sane */
    FCLAW_ASSERT (domain->pp_owned);

    /* check input parameters */
    FCLAW_ASSERT (0 <= blockno && blockno < domain->num_blocks);
    FCLAW_ASSERT (p4est->first_local_tree <= (p4est_topidx_t) blockno);
    FCLAW_ASSERT ((p4est_topidx_t) blockno <= p4est->last_local_tree);
#ifdef FCLAW_ENABLE_DEBUG
    block = domain->blocks + blockno;
#endif
    FCLAW_ASSERT (0 <= patchno && patchno < block->num_patches);

    /* access p4est tree and element for indexing into mesh */
    tree = p4est_tree_array_index (p4est->trees, (p4est_topidx_t) blockno);
    totalleaf = tree->quadrants_offset + (p4est_locidx_t) patchno;
    FCLAW_ASSERT (0 <= totalleaf && totalleaf < p4est->local_num_quadrants);

    /* access face number of the neighbor */
    FCLAW_ASSERT (0 <= faceno && faceno < P4EST_FACES);
    qtfi = (int) mesh->quad_to_face[P4EST_FACES * totalleaf + faceno];
    FCLAW_ASSERT (qtfi >= -num_orient && qtfi < num_orient * (1 + P4EST_HALF));
    FCLAW_ASSERT ((normal_out[faceno] ^
                   normal_out[(qtfi + num_orient) % P4EST_FACES]) ==
                  ((faceno + qtfi) & 1));

    /* we return true if the last bit of the two face numbers differs */
    return (faceno + qtfi) & 1;
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
                                       rproc + 0, rblockno, rpatchno + 0);
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

    P4EST_ASSERT (0 <= *faceno && *faceno < P4EST_FACES);
    P4EST_ASSERT (0 <= *rfaceno && *rfaceno < P4EST_HALF * P4EST_FACES);
    P4EST_ASSERT (0 <= orientation && orientation < P4EST_HALF);

    *faceno = *rfaceno % P4EST_FACES;
    *rfaceno = iface + P4EST_FACES * orientation;

    P4EST_ASSERT (0 <= *faceno && *faceno < P4EST_FACES);
    P4EST_ASSERT (0 <= *rfaceno && *rfaceno < P4EST_HALF * P4EST_FACES);
}

void
fclaw2d_patch_face_transformation (int faceno, int rfaceno, int ftransform[])
{
    p4est_expand_face_transform (faceno, rfaceno, ftransform);
    FCLAW_ASSERT (fclaw2d_patch_face_transformation_valid (ftransform));
}

void
fclaw2d_patch_face_transformation_block (int ftransform[], int sameblock)
{
    FCLAW_ASSERT (fclaw2d_patch_face_transformation_valid (ftransform));
    if (sameblock)
    {
        ftransform[8] |= 4;
    }
    else
    {
        ftransform[8] &= ~4;
    }
    FCLAW_ASSERT (fclaw2d_patch_face_transformation_valid (ftransform));
}

void
fclaw2d_patch_face_transformation_intra (int ftransform[])
{
    memset (ftransform, 0, 9 * sizeof (int));
    ftransform[8] = 4;
}

#ifndef P4_TO_P8
static const int ftransform_max[9] = { 1, 0, 1, 1, 0, 1, 1, 0, 7 };
#else
static const int ftransform_max[9] = { 2, 2, 2, 2, 2, 2, 1, 1, 7 };
#endif

int
fclaw2d_patch_face_transformation_valid (const int ftransform[])
{
    int i;

    for (i = 0; i < 9; ++i) {
        if (ftransform[i] < 0 || ftransform[i] > ftransform_max[i]) {
            return 0;
        }
    }
    return 1;
}

void
fclaw2d_patch_transform_face (fclaw2d_patch_t * ipatch,
                              fclaw2d_patch_t * opatch,
                              const int ftransform[],
                              int mx, int my,
#ifdef P4_TO_P8
                              int mz,
#endif
                              int based, int *i, int *j
#ifdef P4_TO_P8
                            , int *k
#endif
                             )
{
    double iwidth;

    FCLAW_ASSERT (ipatch->level == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && ipatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);
#ifdef P4_TO_P8
    FCLAW_ASSERT (ipatch->zlower >= 0. && ipatch->zlower < 1.);
    FCLAW_ASSERT (opatch->zlower >= 0. && opatch->zlower < 1.);
#endif

    FCLAW_ASSERT (mx >= 1 && my >= 1);
#ifdef P4_TO_P8
    FCLAW_ASSERT (mz >= 1);
#endif
    FCLAW_ASSERT (based == 0 || based == 1);

    FCLAW_ASSERT (fclaw2d_patch_face_transformation_valid (ftransform));

#if 0
#ifndef P4_TO_P8
    printf ("Test I: IP %g %g %d FT %d %d %d %d %d %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level,
            ftransform[0], ftransform[2], ftransform[3], ftransform[5],
            ftransform[6], ftransform[8], mx, *i, *j, based);
#endif
#endif

    /* work with doubles -- exact for integers up to 52 bits of precision */
    iwidth = (double) (1 << ipatch->level);

    if (ftransform[8] & 4)
    {
        /* The two patches are in the same block.  ftransform is not used */
        *i += (int) ((ipatch->xlower - opatch->xlower) * iwidth * (double) mx);
        *j += (int) ((ipatch->ylower - opatch->ylower) * iwidth * (double) my);
#ifdef P4_TO_P8
        *k += (int) ((ipatch->zlower - opatch->zlower) * iwidth * (double) mz);
#endif
    }
    else
    {
        const int *my_axis = &ftransform[0];
        const int *target_axis = &ftransform[3];
        const int *edge_reverse = &ftransform[6];
        double my_xyz[P4EST_DIM], target_xyz[P4EST_DIM];
        double Rmxmymz[P4EST_DIM];
#ifdef FCLAW_ENABLE_DEBUG
        int mxmymz[P4EST_DIM];

        /* make mx, my and mz indexable */
        mxmymz[0] = mx;
        mxmymz[1] = my;
#ifdef P4_TO_P8
        mxmymz[2] = mz;
#endif
#endif

        /* make gridsize indexable */
        Rmxmymz[0] = iwidth * (double) mx;
        Rmxmymz[1] = iwidth * (double) my;
#ifdef P4_TO_P8
        Rmxmymz[2] = iwidth * (double) mz;
#endif

        /* the reference cube is stretched to mx times my units */
        my_xyz[0] = ipatch->xlower * Rmxmymz[0] + *i - .5 * based;
        my_xyz[1] = ipatch->ylower * Rmxmymz[1] + *j - .5 * based;
#ifdef P4_TO_P8
        my_xyz[2] = ipatch->zlower * Rmxmymz[2] + *k - .5 * based;
#endif

        /* transform transversal directions */
        FCLAW_ASSERT (mxmymz[my_axis[0]] == mxmymz[target_axis[0]]);
        target_xyz[target_axis[0]] =
            !edge_reverse[0] ? my_xyz[my_axis[0]] : Rmxmymz[my_axis[0]] -
            my_xyz[my_axis[0]];
#ifdef P4_TO_P8
        FCLAW_ASSERT (mxmymz[my_axis[1]] == mxmymz[target_axis[1]]);
        target_xyz[target_axis[1]] =
            !edge_reverse[1] ? my_xyz[my_axis[1]] : Rmxmymz[my_axis[1]] -
            my_xyz[my_axis[1]];
#endif

        /* transform normal direction */
        FCLAW_ASSERT (mxmymz[my_axis[2]] == mxmymz[target_axis[2]]);
        switch (edge_reverse[2])
        {
        case 0:
            target_xyz[target_axis[2]] = -my_xyz[my_axis[2]];
            break;
        case 1:
            target_xyz[target_axis[2]] =
                my_xyz[my_axis[2]] + Rmxmymz[my_axis[2]];
            break;
        case 2:
            target_xyz[target_axis[2]] =
                my_xyz[my_axis[2]] - Rmxmymz[my_axis[2]];
            break;
        case 3:
            target_xyz[target_axis[2]] =
                2. * Rmxmymz[my_axis[2]] - my_xyz[my_axis[2]];
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* transform back to integer coordinates: this is exact */
        *i = (int) (target_xyz[0] - opatch->xlower * Rmxmymz[0] + .5 * based);
        *j = (int) (target_xyz[1] - opatch->ylower * Rmxmymz[1] + .5 * based);
#ifdef P4_TO_P8
        *k = (int) (target_xyz[2] - opatch->zlower * Rmxmymz[2] + .5 * based);
#endif
    }

#if 0
#ifndef P4_TO_P8
    printf ("Test O: IP %g %g IJ %d %d\n",
            opatch->xlower, opatch->ylower, *i, *j);
#endif
#endif
}

void
fclaw2d_patch_transform_face2 (fclaw2d_patch_t * ipatch,
                               fclaw2d_patch_t * opatch,
                               const int ftransform[],
                               int mx, int my,
#ifdef P4_TO_P8
                               int mz,
#endif
                               int based, int i[], int j[]
#ifdef P4_TO_P8
                             , int k[]
#endif
                              )
{
    int di, dj;
#ifdef P4_TO_P8
    int dk;
#endif
    double owidth;

    FCLAW_ASSERT (ipatch->level + 1 == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && opatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);
#ifdef P4_TO_P8
    FCLAW_ASSERT (ipatch->zlower >= 0. && ipatch->zlower < 1.);
    FCLAW_ASSERT (opatch->zlower >= 0. && opatch->zlower < 1.);
#endif

    FCLAW_ASSERT (mx >= 1 && my >= 1);
#ifdef P4_TO_P8
    FCLAW_ASSERT (mz >= 1);
#endif
    FCLAW_ASSERT (based == 0 || based == 1);

    FCLAW_ASSERT (fclaw2d_patch_face_transformation_valid (ftransform));

#if 0
#ifndef P4_TO_P8
    printf ("Test I: IP %g %g %d FT %d %d %d %d %d %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level,
            ftransform[0], ftransform[2], ftransform[3], ftransform[5],
            ftransform[6], ftransform[8], mx, *i, *j, based);
#endif
#endif

    /* work with doubles -- exact for integers up to 52 bits of precision */
    owidth = (double) (1 << opatch->level);

    if (ftransform[8] & 4)
    {
        int kx, ky, kz;

        /* The two patches are in the same block.  ftransform is undefined */
        di = based + (int)
            ((ipatch->xlower - opatch->xlower) * owidth * (double) mx +
             2. * (*i - based));
        dj = based +
            (int) ((ipatch->ylower - opatch->ylower) * owidth * (double) my +
                   2. * (*j - based));
#ifdef P4_TO_P8
        dk = based + (int)
            ((ipatch->zlower - opatch->zlower) * owidth * (double) mz +
             2. * (*k - based));
#else
        kz = 0;
#endif

        /* In the same block, the order of child cells is canonical */
#ifdef P4_TO_P8
        for (kz = 0; kz < 2; ++kz) {
#if 0
        }
#endif
#endif
        for (ky = 0; ky < 2; ++ky)
        {
            for (kx = 0; kx < 2; ++kx)
            {
                i[4 * kz + 2 * ky + kx] = di + kx;
                j[4 * kz + 2 * ky + kx] = dj + ky;
#ifdef P4_TO_P8
                k[4 * kz + 2 * ky + kx] = dk + kz;
#endif
            }
        }
#ifdef P4_TO_P8
#if 0
        {
#endif
        }
#endif
    }
    else
    {
        int l, is[3], ik[3], ib[3], in[3];
        const int *my_axis = &ftransform[0];
        const int *target_axis = &ftransform[3];
        const int *edge_reverse = &ftransform[6];
        double my_xyz[P4EST_DIM], target_xyz[P4EST_DIM];
        double Rmxmymz[P4EST_DIM];
#ifdef FCLAW_ENABLE_DEBUG
        int mxmymz[P4EST_DIM];

        /* make mx, my and mz indexable */
        mxmymz[0] = mx;
        mxmymz[1] = my;
#ifdef P4_TO_P8
        mxmymz[2] = mz;
#endif
#endif
        /* make gridsize indexable */
        Rmxmymz[0] = owidth * (double) mx;
        Rmxmymz[1] = owidth * (double) my;
#ifdef P4_TO_P8
        Rmxmymz[2] = owidth * (double) mz;
#endif

        /* the reference cube is stretched to mx times my units */
        my_xyz[0] = ipatch->xlower * Rmxmymz[0] + 2. * (*i + .5 - based);
        my_xyz[1] = ipatch->ylower * Rmxmymz[1] + 2. * (*j + .5 - based);
#ifdef P4_TO_P8
        my_xyz[2] = ipatch->zlower * Rmxmymz[2] + 2. * (*k + .5 - based);
#else
        is[2] = ik[1] = ib[1] = in[2] = 0;
#endif

        /* transform transversal directions */
        FCLAW_ASSERT (mxmymz[target_axis[0]] == mxmymz[my_axis[0]]);
        target_xyz[target_axis[0]] =
            !edge_reverse[0] ? my_xyz[my_axis[0]] : Rmxmymz[my_axis[0]] -
            my_xyz[my_axis[0]];
#ifdef P4_TO_P8
        FCLAW_ASSERT (mxmymz[target_axis[1]] == mxmymz[my_axis[1]]);
        target_xyz[target_axis[1]] =
            !edge_reverse[1] ? my_xyz[my_axis[1]] : Rmxmymz[my_axis[1]] -
            my_xyz[my_axis[1]];
#endif

        /* transform normal direction */
        switch (edge_reverse[2])
        {
        case 0:
            target_xyz[target_axis[2]] = -my_xyz[my_axis[2]];
            break;
        case 1:
            target_xyz[target_axis[2]] =
                my_xyz[my_axis[2]] + Rmxmymz[my_axis[2]];
            break;
        case 2:
            target_xyz[target_axis[2]] =
                my_xyz[my_axis[2]] - Rmxmymz[my_axis[2]];
            break;
        case 3:
            target_xyz[target_axis[2]] =
                2. * Rmxmymz[my_axis[2]] - my_xyz[my_axis[2]];
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }

        /* move back into integer coordinates: this is exact */
        di = (int) (target_xyz[0] - opatch->xlower * Rmxmymz[0]) + based - 1;
        dj = (int) (target_xyz[1] - opatch->ylower * Rmxmymz[1]) + based - 1;
#ifdef P4_TO_P8
        dk = (int) (target_xyz[2] - opatch->zlower * Rmxmymz[2]) + based - 1;
#endif

        /* Run through the child cells in order of the small patch */
#ifdef P4_TO_P8
        for (is[2] = 0; is[2] < 2; ++is[2]) {
#if 0
        }
#endif
#endif
        for (is[1] = 0; is[1] < 2; ++is[1])
        {
            for (is[0] = 0; is[0] < 2; ++is[0])
            {
                /* Compute small patch coordinate in (transverse, normal) order */
                for (l = 0; l < 3; ++l)
                {
#if (P4EST_DIM == 2)
                    if (l == 1) continue;
#endif
                    ik[l] = is[target_axis[l]];
                }

                /* Transform into big patch (transverse, normal) coordinate */
                ib[0] = !edge_reverse[0] ? ik[0] : !ik[0];
#ifdef P4_TO_P8
                ib[1] = !edge_reverse[1] ? ik[1] : !ik[1];
#endif
                ib[2] = (edge_reverse[2] == 1
                      || edge_reverse[2] == 2) ? ik[2] : !ik[2];

                /* Compute coordinate relative to the big patch */
                for (l = 0; l < 3; ++l)
                {
#if (P4EST_DIM == 2)
                    if (l == 1) continue;
#endif
                    in[my_axis[l]] = ib[l];
                }

                /* assign value in proper place */
                i[4 * in[2] + 2 * in[1] + in[0]] = di + is[0];
                j[4 * in[2] + 2 * in[1] + in[0]] = dj + is[1];
#ifdef P4_TO_P8
                k[4 * in[2] + 2 * in[1] + in[0]] = dk + is[2];
#endif
            }
        }
#ifdef P4_TO_P8
#if 0
        {
#endif
        }
#endif
    }

#if 0
#ifndef P4_TO_P8
    printf ("Test O: OP %g %g %d I %d %d %d %d J %d %d %d %d\n",
            opatch->xlower, opatch->ylower, opatch->level,
            i[0], i[1], i[2], i[3], j[0], j[1], j[2], j[3]);
#endif
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
     * than four (2D) or eight (3D) blocks meet at a corner */
    if (qid >= 0)
    {
        FCLAW_ASSERT (0 <= qid);
        if (qid >= mesh->local_num_quadrants + mesh->ghost_num_quadrants)
        {
            /* This is an inter-tree (face, edge or corner) corner neighbor */
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
                /* At least a five/nine-corner, which is currently not supported */
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
        *rblockno = blockno;
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
                                int mx, int my,
#ifdef P4_TO_P8
                                int mz,
#endif
                                int based, int *i, int *j
#ifdef P4_TO_P8
                              , int *k
#endif
                               )
{
    double Rmxmymz[P4EST_DIM], xshift, yshift;
#ifdef P4_TO_P8
    double zshift;
#endif

    FCLAW_ASSERT (ipatch->level == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && ipatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);
#ifdef P4_TO_P8
    FCLAW_ASSERT (ipatch->zlower >= 0. && ipatch->zlower < 1.);
    FCLAW_ASSERT (opatch->zlower >= 0. && opatch->zlower < 1.);
#endif

    FCLAW_ASSERT (mx >= 1 && my >= 1);
#ifdef P4_TO_P8
    FCLAW_ASSERT (mz >= 1);
#endif
    FCLAW_ASSERT (based == 0 || based == 1);

#if 0
#ifndef P4_TO_P8
    printf ("Test I: IP %g %g %d IC %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level, icorner,
            mx, *i, *j, based);
#endif
#endif

    /* Work with doubles -- exact for integers up to 52 bits of precision */
    Rmxmymz[0] = (double) (1 << ipatch->level) * (double) mx;
    Rmxmymz[1] = (double) (1 << ipatch->level) * (double) my;
#ifdef P4_TO_P8
    Rmxmymz[2] = (double) (1 << ipatch->level) * (double) mz;
#endif
    if (!is_block_boundary)
    {
        /* The lower left coordinates are with respect to the same origin */
        xshift = yshift = 0.;
#ifdef P4_TO_P8
        zshift = 0.;
#endif
    }
    else
    {
        /* We need to add/substract the shift due to translation of one block */
        if ((icorner & 1) == 0)
        {
            /* This corner is on the left face of the patch */
            /* verify the blocks are corner-neighbors with the next assertions */
            FCLAW_ASSERT (fabs (opatch->xupper - 1.) < SC_1000_EPS);
            xshift = +1.;
        }
        else
        {
            /* This corner is on the right face of the patch */
            FCLAW_ASSERT (fabs (opatch->xlower) < SC_1000_EPS);
            xshift = -1.;
        }
        if ((icorner & 2) == 0)
        {
            /* This corner is on the front face of the patch */
            FCLAW_ASSERT (fabs (opatch->yupper - 1.) < SC_1000_EPS);
            yshift = +1.;
        }
        else
        {
            /* This corner is on the back face of the patch */
            FCLAW_ASSERT (fabs (opatch->ylower) < SC_1000_EPS);
            yshift = -1.;
        }
#ifdef P4_TO_P8
        if ((icorner & 4) == 0)
        {
            /* This corner is on the bottom face of the patch */
            FCLAW_ASSERT (fabs (opatch->zupper - 1.) < SC_1000_EPS);
            zshift = +1.;
        }
        else
        {
            /* This corner is on the top face of the patch */
            FCLAW_ASSERT (fabs (opatch->zlower) < SC_1000_EPS);
            zshift = -1.;
        }
#endif
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    *i += (int) ((ipatch->xlower - opatch->xlower + xshift) * Rmxmymz[0]);
    *j += (int) ((ipatch->ylower - opatch->ylower + yshift) * Rmxmymz[1]);
#ifdef P4_TO_P8
    *k += (int) ((ipatch->zlower - opatch->zlower + zshift) * Rmxmymz[2]);
#endif
}

void
fclaw2d_patch_transform_corner2 (fclaw2d_patch_t * ipatch,
                                 fclaw2d_patch_t * opatch,
                                 int icorner, int is_block_boundary,
                                 int mx, int my,
#ifdef P4_TO_P8
                                 int mz,
#endif
                                 int based, int i[], int j[]
#ifdef P4_TO_P8
                               , int k[]
#endif
                                )
{
    int kt, kn, ks;
    int di, dj;
#ifdef P4_TO_P8
    int dk;
    double zshift;
#endif
    double Rmxmymz[P4EST_DIM], xshift, yshift;

    FCLAW_ASSERT (ipatch->level + 1 == opatch->level);
    FCLAW_ASSERT (0 <= ipatch->level && opatch->level < P4EST_MAXLEVEL);
    FCLAW_ASSERT (ipatch->xlower >= 0. && ipatch->xlower < 1.);
    FCLAW_ASSERT (opatch->xlower >= 0. && opatch->xlower < 1.);
    FCLAW_ASSERT (ipatch->ylower >= 0. && ipatch->ylower < 1.);
    FCLAW_ASSERT (opatch->ylower >= 0. && opatch->ylower < 1.);
#ifdef P4_TO_P8
    FCLAW_ASSERT (ipatch->zlower >= 0. && ipatch->zlower < 1.);
    FCLAW_ASSERT (opatch->zlower >= 0. && opatch->zlower < 1.);
#endif

    FCLAW_ASSERT (mx >= 1 && my >= 1);
#ifdef P4_TO_P8
    FCLAW_ASSERT (mz >= 1);
#endif
    FCLAW_ASSERT (based == 0 || based == 1);

#if 0
#ifndef P4_TO_P8
    printf ("Test I: IP %g %g %d IC %d MX %d IJ %d %d BS %d\n",
            ipatch->xlower, ipatch->ylower, ipatch->level, icorner,
            mx, *i, *j, based);
#endif
#endif

    /* work with doubles -- exact for integers up to 52 bits of precision */
    Rmxmymz[0] = (double) (1 << opatch->level) * (double) mx;
    Rmxmymz[1] = (double) (1 << opatch->level) * (double) my;
#ifdef P4_TO_P8
    Rmxmymz[2] = (double) (1 << opatch->level) * (double) mz;
#endif
    if (!is_block_boundary)
    {
        /* The lower left coordinates are with respect to the same origin */
        xshift = yshift = 0.;
#ifdef P4_TO_P8
        zshift = 0.;
#endif
    }
    else
    {
        /* We need to add/substract the shift due to translation of one block */
        if ((icorner & 1) == 0)
        {
            /* This corner is on the left face of the patch */
            /* verify the blocks are corner-neighbors with the next assertions */
            FCLAW_ASSERT (fabs (opatch->xupper - 1.) < SC_1000_EPS);
            xshift = +1.;
        }
        else
        {
            /* This corner is on the right face of the patch */
            FCLAW_ASSERT (fabs (opatch->xlower) < SC_1000_EPS);
            xshift = -1.;
        }
        if ((icorner & 2) == 0)
        {
            /* This corner is on the front face of the patch */
            FCLAW_ASSERT (fabs (opatch->yupper - 1.) < SC_1000_EPS);
            yshift = +1.;
        }
        else
        {
            /* This corner is on the back face of the patch */
            FCLAW_ASSERT (fabs (opatch->ylower) < SC_1000_EPS);
            yshift = -1.;
        }
#ifdef P4_TO_P8
        if ((icorner & 4) == 0)
        {
            /* This corner is on the bottom face of the patch */
            FCLAW_ASSERT (fabs (opatch->zupper - 1.) < SC_1000_EPS);
            zshift = +1.;
        }
        else
        {
            /* This corner is on the top face of the patch */
            FCLAW_ASSERT (fabs (opatch->zlower) < SC_1000_EPS);
            zshift = -1.;
        }
#endif
    }

    /* The two patches are in the same block, or in a different block
     * that has a coordinate system with the same orientation */
    di = based
        + (int) ((ipatch->xlower - opatch->xlower + xshift) * Rmxmymz[0] +
                 2. * (*i - based));
    dj = based
        + (int) ((ipatch->ylower - opatch->ylower + yshift) * Rmxmymz[1] +
                 2. * (*j - based));
#ifdef P4_TO_P8
    dk = based
        + (int) ((ipatch->zlower - opatch->zlower + zshift) * Rmxmymz[2] +
                 2. * (*k - based));
#else
    ks = 0;
#endif

    /* Without any rotation, the order of child cells is canonical */
#ifdef P4_TO_P8
    for (ks = 0; ks < 2; ++ks) {
#if 0
    }
#endif
#endif
    for (kt = 0; kt < 2; ++kt)
    {
        for (kn = 0; kn < 2; ++kn)
        {
            i[4 * ks + 2 * kt + kn] = di + kn;
            j[4 * ks + 2 * kt + kn] = dj + kt;
#ifdef P4_TO_P8
            k[4 * ks + 2 * kt + kn] = dk + ks;
#endif
        }
    }
#ifdef P4_TO_P8
#if 0
    {
#endif
    }
#endif
}

void
fclaw2d_domain_set_refinement (fclaw2d_domain_t * domain,
                               int smooth_refine, int smooth_level,
                               int coarsen_delay)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    FCLAW_ASSERT (coarsen_delay >= 0);

    domain->p.smooth_refine = smooth_refine;
    domain->p.smooth_level = smooth_level;
    p4est_wrap_set_coarsen_delay (wrap, coarsen_delay, 0);
}

void
fclaw2d_domain_set_partitioning (fclaw2d_domain_t * domain,
                                 int partition_for_coarsening)
{
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    p4est_wrap_set_partitioning (wrap, partition_for_coarsening);
}

void
fclaw2d_patch_mark_refine (fclaw2d_domain_t * domain, int blockno,
                           int patchno)
{
    fclaw2d_patch_t *patch;
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;

    /* bump target level up */
    patch = fclaw2d_domain_get_patch (domain, blockno, patchno);
    /* Only tag for refinement if current level is greater
       or equal to than smooth level.  Any refinement that
       happens below smooth_level is done only for
       balancing purposes */

#if 0
    patch->target_level =
        patch->level >= domain->p.smooth_level ?
        patch->level + 1 : patch->level;
#endif

    patch->target_level = patch->level + 1;
    patch->target_level = SC_MIN (patch->target_level, P4EST_QMAXLEVEL);

    /* if we do smooth refinement, all marking is done inside adapt */
    if (!domain->p.smooth_refine)
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
    if (!domain->p.smooth_refine)
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

struct fclaw2d_domain_indirect
{
    int ready;
    fclaw2d_domain_t *domain;
    fclaw2d_domain_exchange_t *e;
};

static void
indirect_encode (p4est_ghost_t * ghost, int mpirank,
                 int *rproc, int *rpatchno)
{
    p4est_quadrant_t *g;

    FCLAW_ASSERT (0 <= *rproc && *rproc < ghost->mpisize);
    if (*rproc == mpirank)
    {
        /* includes the case FCLAW2D_PATCH_BOUNDARY */
        *rproc = *rpatchno = -1;
    }
    else
    {
        FCLAW_ASSERT (ghost->proc_offsets[*rproc] <= *rpatchno &&
                      *rpatchno < ghost->proc_offsets[*rproc + 1]);
        g = p4est_quadrant_array_index (&ghost->ghosts, *rpatchno);
        *rpatchno = (int) g->p.piggy3.local_num;
        FCLAW_ASSERT (*rpatchno >= 0);
    }
}

static void
indirect_match_face (int *pi,
                     int **rproc, int **rblockno, int **rpatchno,
                     int **rfaceno)
{
    FCLAW_ASSERT (pi != NULL);
    FCLAW_ASSERT (rproc != NULL && rblockno != NULL);
    FCLAW_ASSERT (rpatchno != NULL && rfaceno != NULL);

    *rproc = pi;
    *rblockno = pi + P4EST_HALF;
    *rpatchno = pi + P4EST_HALF + 1;
    *rfaceno = pi + P4EST_HALF + 1 + P4EST_HALF;
}

#ifdef P4_TO_P8
static void
indirect_match_edge (int *pi,
                     int **rproc, int **rblockno, int **rpatchno,
                     int **rfaceno)
{
    FCLAW_ASSERT (pi != NULL);
    FCLAW_ASSERT (rproc != NULL && rblockno != NULL);
    FCLAW_ASSERT (rpatchno != NULL && rfaceno != NULL);

    *rproc = pi;
    *rblockno = pi + 2;
    *rpatchno = pi + 2 + 1;
    *rfaceno = pi + 2 + 1 + 2;
}
#endif

static void
indirect_match_corner (int *pi,
                       int **rproc, int **rblockno, int **rpatchno,
                       int **rcornerno)
{
    FCLAW_ASSERT (pi != NULL);
    FCLAW_ASSERT (rproc != NULL && rblockno != NULL);
    FCLAW_ASSERT (rpatchno != NULL && rcornerno != NULL);

    *rproc = pi;
    *rblockno = pi + 1;
    *rpatchno = pi + 2;
    *rcornerno = pi + 3;
}

fclaw2d_domain_indirect_t *
fclaw2d_domain_indirect_begin (fclaw2d_domain_t * domain)
{
    int num_exc;
    int neall, nb, ne, np;
    int face, corner;
#ifdef P4_TO_P8
    int edge;
#endif
    int *pbdata, *pi;
    int *rproc, *rblockno, *rpatchno, *rfaceno, *rcornerno;
#ifdef P4_TO_P8
    int *redgeno;
#endif
    int has_corner;
    int face_info_size, corner_info_size, info_size;
#ifdef P4_TO_P8
    int edge_info_size;
#endif
    size_t data_size;
    int sf;
    fclaw2d_block_t *block;
    fclaw2d_domain_indirect_t *ind;
    fclaw2d_patch_relation_t prel;
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_ghost_t *ghost = p4est_wrap_get_ghost (wrap);

    num_exc = domain->num_exchange_patches;
    face_info_size = 2 + 2 * P4EST_HALF;
#ifdef P4_TO_P8
    edge_info_size = 2 + 2 * 2;
#endif
    corner_info_size = 4;
    info_size =
        P4EST_FACES * face_info_size + P4EST_CHILDREN * corner_info_size;
#ifdef P4_TO_P8
    info_size += P8EST_EDGES * edge_info_size;
#endif
    data_size = info_size * sizeof (int);

    /* allocate internal state for this operation */
    ind = FCLAW_ALLOC_ZERO (fclaw2d_domain_indirect_t, 1);
    ind->domain = domain;
    ind->e = fclaw2d_domain_allocate_before_exchange (domain, data_size);

    /* loop through exchange patches and fill their neighbor information */
    pbdata = pi = FCLAW_ALLOC (int, num_exc * info_size);
    for (neall = 0, nb = 0; nb < domain->num_blocks; ++nb)
    {
        block = domain->blocks + nb;
        for (ne = 0; ne < block->num_exchange_patches; ++ne, ++neall)
        {
            ind->e->patch_data[neall] = (void *) pi;
            np = (int) (block->exchange_patches[ne] - block->patches);
            FCLAW_ASSERT (0 <= np && np < block->num_patches);
            FCLAW_ASSERT (block->exchange_patches[ne] ==
                          domain->exchange_patches[neall]);
            for (face = 0; face < P4EST_FACES; ++face)
            {
                indirect_match_face (pi, &rproc, &rblockno, &rpatchno,
                                     &rfaceno);
                prel =
                    fclaw2d_patch_face_neighbors (domain, nb, np, face, rproc,
                                                  rblockno, rpatchno,
                                                  rfaceno);

                /* obtain proper ghost patch numbers for the receiver */
                indirect_encode (ghost, domain->mpirank,
                                 &rproc[0], &rpatchno[0]);
                if (prel == FCLAW2D_PATCH_HALFSIZE)
                {
                    for (sf = 1; sf < P4EST_HALF; sf++)
                    {
                        indirect_encode
                            (ghost, domain->mpirank, &rproc[sf],
                             &rpatchno[sf]);
                    }
                    *rfaceno |= 1 << 26;
                }
                else if (prel == FCLAW2D_PATCH_DOUBLESIZE)
                {
                    *rfaceno |= 1 << 27;
                }
                pi += face_info_size;
            }

            for (corner = 0; corner < P4EST_CHILDREN; ++corner)
            {
                indirect_match_corner (pi, &rproc, &rblockno, &rpatchno,
                                       &rcornerno);
                has_corner =
                    fclaw2d_patch_corner_neighbors (domain, nb, np, corner,
                                                    rproc, rblockno, rpatchno,
                                                    rcornerno, &prel);

                if (has_corner)
                {
                    /* obtain proper ghost patch numbers for the receiver */
                    indirect_encode (ghost, domain->mpirank,
                                     &rproc[0], &rpatchno[0]);
                    if (prel == FCLAW2D_PATCH_HALFSIZE)
                    {
                        *rcornerno |= 1 << 26;
                    }
                    else if (prel == FCLAW2D_PATCH_DOUBLESIZE)
                    {
                        *rcornerno |= 1 << 27;
                    }
                }
                else
                {
                    rproc[0] = rpatchno[0] = -1;
                    *rcornerno = corner;        /* rcornerno == -1 messes with & operator */
                }
                pi += corner_info_size;
            }

#ifdef P4_TO_P8
            for(edge = 0; edge < P8EST_EDGES; ++edge)
            {
                indirect_match_edge (pi, &rproc, &rblockno, &rpatchno,
                                     &redgeno);

                fclaw3d_patch_edge_neighbors (domain, nb, np, edge, rproc,
                                              rblockno, rpatchno,
                                              redgeno, &prel);
                FCLAW_ASSERT (0 <= *rblockno && *rblockno < domain->num_blocks);

                if(prel != FCLAW2D_PATCH_BOUNDARY)
                {
                    /* obtain proper ghost patch numbers for the receiver */
                    indirect_encode (ghost, domain->mpirank,
                                     &rproc[0], &rpatchno[0]);
                    if (prel == FCLAW2D_PATCH_HALFSIZE)
                    {
                        indirect_encode
                            (ghost, domain->mpirank, &rproc[1],
                             &rpatchno[1]);
                        *redgeno |= 1 << 26;
                    }
                    else if (prel == FCLAW2D_PATCH_DOUBLESIZE)
                    {
                        *redgeno |= 1 << 27;
                    }
                }
                else
                {
                    rproc[0] = rpatchno[0] = -1;
                    rproc[1] = rpatchno[1] = -1;
                    *redgeno = edge;        /* rcornerno == -1 messes with & operator */
                }
                pi += edge_info_size;
            }
#endif 
        }
    }
    FCLAW_ASSERT ((int) (pi - pbdata) == neall * info_size);
    FCLAW_ASSERT (neall == num_exc);

    /* post messages */
    fclaw2d_domain_ghost_exchange_begin (domain, ind->e,
                                         0, domain->global_maxlevel);

    /* it is safe to destroy the sent data now; they have been copied */
    FCLAW_FREE (pbdata);

    return ind;
}

/* These static functions are unused in 3D as long as the 3D case is not
 * implemented.
 */
static uint64_t
pli_make_key (int p, int rpatchno)
{
    FCLAW_ASSERT (p >= 0 && rpatchno >= 0);
    return (((uint64_t) p) << 32) + rpatchno;
}

static unsigned
pli_hash_function (const void *v, const void *u)
{
    const uint64_t ui1 = *(uint64_t *) v;
    uint32_t i1, i2, i3;

    i1 = (uint32_t) ui1;
    i2 = (uint32_t) (ui1 >> 16);
    i3 = (uint32_t) (ui1 >> 32);
    sc_hash_mix (i1, i2, i3);
    sc_hash_final (i1, i2, i3);

    return (unsigned) i1;
}

static int
pli_equal_function (const void *v1, const void *v2, const void *u)
{
    const uint64_t ui1 = *(uint64_t *) v1;
    const uint64_t ui2 = *(uint64_t *) v2;

    return ui1 == ui2;
}

static int
pli_make_ng (void **found, uint64_t * pli_keys)
{
    return (int) ((uint64_t *) * found - pli_keys);
}

static int
indirect_decode (sc_hash_t * pli_hash, uint64_t * pli_keys,
                 int mpisize, int mpirank, int *rproc, int *rpatchno)
{
    int good = 0;
    int retval;
    uint64_t akey;
    void **found;

    if (*rproc == mpirank)
    {
        *rproc = *rpatchno = -1;
    }
    if (*rproc != -1)
    {
        FCLAW_ASSERT (0 <= *rproc && *rproc < mpisize);
        akey = pli_make_key (*rproc, *rpatchno);
        retval = sc_hash_lookup (pli_hash, &akey, &found);
        if (retval)
        {
            *rpatchno = pli_make_ng (found, pli_keys);
            good = 1;
        }
        else
        {
            *rproc = *rpatchno = -1;
        }
    }
    else
    {
        FCLAW_ASSERT (*rpatchno == -1);
    }

    return good;
}

void
fclaw2d_domain_indirect_end (fclaw2d_domain_t * domain,
                             fclaw2d_domain_indirect_t * ind)
{
    int ndgp;
    int good, good2;
    int p, ng;
#ifdef FCLAW_ENABLE_DEBUG
    int gprev;
#endif
    int gpatch;
    int face, corner;
    int *rproc, *rblockno, *rpatchno, *rfaceno, *rcornerno;
    int *pi;
    int face_info_size, corner_info_size;
    int sf;
#ifdef P4_TO_P8
    int edge;
    int *redgeno;
    int edge_info_size;
#endif
    uint64_t *pli_keys, *plik;
    sc_hash_t *pli_hash;
    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_ghost_t *ghost = p4est_wrap_get_ghost (wrap);
    p4est_quadrant_t *g;

    FCLAW_ASSERT (ind != NULL && !ind->ready);
    FCLAW_ASSERT (domain == ind->domain);
    ndgp = domain->num_ghost_patches;

    /* prepare lookup logic for processor-local ghost indices */
    pli_keys = plik = FCLAW_ALLOC (uint64_t, ndgp);
    pli_hash = sc_hash_new (pli_hash_function, pli_equal_function,
                            NULL, NULL);
    for (ng = 0, p = 0; p < domain->mpisize; ++p)
    {
#ifdef FCLAW_ENABLE_DEBUG
        gprev = -1;
#endif
        for (; ng < (int) ghost->proc_offsets[p + 1]; ++ng)
        {
            g = p4est_quadrant_array_index (&ghost->ghosts, ng);
            gpatch = (int) g->p.piggy3.local_num;
            FCLAW_ASSERT (gprev < gpatch);
            *plik = pli_make_key (p, gpatch);
            SC_EXECUTE_ASSERT_TRUE
                (sc_hash_insert_unique (pli_hash, plik, NULL));
#ifdef FCLAW_ENABLE_DEBUG
            gprev = gpatch;
#endif
            ++plik;
        }
    }
    FCLAW_ASSERT ((int) (plik - pli_keys) == ndgp);

    /* receive messages */
    fclaw2d_domain_ghost_exchange_end (domain, ind->e);

    /* go through ghosts a second time, now working on received data */
    face_info_size = 2 + 2 * P4EST_HALF;
    corner_info_size = 4;
#ifdef P4_TO_P8
    edge_info_size = 2 + 2 * 2;
#endif
    for (ng = 0, p = 0; p < domain->mpisize; ++p)
    {
        for (; ng < (int) ghost->proc_offsets[p + 1]; ++ng)
        {
            g = p4est_quadrant_array_index (&ghost->ghosts, ng);
            gpatch = (int) g->p.piggy3.local_num;
            pi = (int *) ind->e->ghost_data[ng];

            /* go through face neighbor patches of this ghost */
            for (face = 0; face < P4EST_FACES; ++face)
            {
                /* access values that were shipped with the ghost */
                indirect_match_face (pi, &rproc, &rblockno, &rpatchno,
                                     &rfaceno);

                /* check first face neighbor */
                FCLAW_ASSERT (rproc[0] != p);
                good = indirect_decode (pli_hash, pli_keys,
                                        domain->mpisize, domain->mpirank,
                                        &rproc[0], &rpatchno[0]);
                FCLAW_ASSERT ((rproc[0] == -1 && rpatchno[0] == -1) ||
                              (0 <= rpatchno[0] && rpatchno[0] < ndgp));
                if (*rfaceno & (1 << 26))
                {
                    /* check remaining of P4EST_HALF halfsize neighbors */
                    for (sf = 1; sf < P4EST_HALF; sf++)
                    {
                        FCLAW_ASSERT (rproc[sf] != p);
                        good2 = indirect_decode (pli_hash, pli_keys,
                                                 domain->mpisize,
                                                 domain->mpirank, &rproc[sf],
                                                 &rpatchno[sf]);
                        good = good || good2;
                        FCLAW_ASSERT ((rproc[sf] == -1 && rpatchno[sf] == -1)
                                      || (0 <= rpatchno[sf]
                                          && rpatchno[sf] < ndgp));
                    }
                }

                /* no match on this face; we pretend a boundary situation */
                if (!good)
                {
                    FCLAW_ASSERT (rproc[0] == -1);
                    FCLAW_ASSERT (rpatchno[0] == -1);
                    for (sf = 1; sf < P4EST_HALF; sf++)
                    {
                        rproc[sf] = rpatchno[sf] = -1;
                    }
                    *rfaceno &= ~(3 << 26);
                }

                /* move to the next face data item */
                pi += face_info_size;
            }

            /* go through corner neighbor patches of this ghost */
            for (corner = 0; corner < P4EST_CHILDREN; ++corner)
            {
                /* access values that were shipped with the ghost */
                indirect_match_corner (pi, &rproc, &rblockno, &rpatchno,
                                       &rcornerno);

                FCLAW_ASSERT (rproc[0] != p);
                good = indirect_decode (pli_hash, pli_keys,
                                        domain->mpisize, domain->mpirank,
                                        &rproc[0], &rpatchno[0]);
                FCLAW_ASSERT ((rproc[0] == -1 && rpatchno[0] == -1) ||
                              (0 <= rpatchno[0] && rpatchno[0] < ndgp));

                /* no match on this corner; we pretend a boundary situation */
                if (!good)
                {
                    FCLAW_ASSERT (rproc[0] == -1);
                    FCLAW_ASSERT (rpatchno[0] == -1);
                    *rcornerno &= ~(3 << 26);
                }

                /* move to the next face data item */
                pi += corner_info_size;
            }

#ifdef P4_TO_P8
            /* go through edge neighbor patches of this ghost */
            for (edge = 0; edge < P8EST_EDGES; ++edge)
            {
                /* access values that were shipped with the ghost */
                indirect_match_edge (pi, &rproc, &rblockno, &rpatchno,
                                     &redgeno);

                /* check first edge neighbor */
                FCLAW_ASSERT (rproc[0] != p);
                good = indirect_decode (pli_hash, pli_keys,
                                        domain->mpisize, domain->mpirank,
                                        &rproc[0], &rpatchno[0]);
                FCLAW_ASSERT ((rproc[0] == -1 && rpatchno[0] == -1) ||
                              (0 <= rpatchno[0] && rpatchno[0] < ndgp));
                if (*redgeno & (1 << 26))
                {
                    /* check other halfsize neighbor */
                    FCLAW_ASSERT (rproc[1] != p);
                    good2 = indirect_decode (pli_hash, pli_keys,
                                             domain->mpisize,
                                             domain->mpirank, &rproc[1],
                                             &rpatchno[1]);
                    good = good || good2;
                    FCLAW_ASSERT ((rproc[1] == -1 && rpatchno[1] == -1)
                                  || (0 <= rpatchno[1]
                                      && rpatchno[1] < ndgp));
                }

                /* no match on this edge; we pretend a boundary situation */
                if (!good)
                {
                    FCLAW_ASSERT (rproc[0] == -1);
                    FCLAW_ASSERT (rpatchno[0] == -1);
                    rproc[1] = rpatchno[1] = -1;
                    *redgeno &= ~(3 << 26);
                }

                /* move to the next edge data item */
                pi += edge_info_size;
            }
#endif
        }
    }
    FCLAW_ASSERT (ng == domain->num_ghost_patches);

    /* hash data is local to this function */
    FCLAW_FREE (pli_keys);
    sc_hash_destroy (pli_hash);

    /* now we allow queries on the ghost data */
    ind->ready = 1;
}

fclaw2d_patch_relation_t
fclaw2d_domain_indirect_face_neighbors (fclaw2d_domain_t * domain,
                                        fclaw2d_domain_indirect_t * ind,
                                        int ghostno, int faceno,
                                        int rproc[P4EST_HALF], int *rblockno,
                                        int rpatchno[P4EST_HALF],
                                        int *rfaceno)
{
    int *pi;
    int *grproc, *grblockno, *grpatchno, *grfaceno;
    int sf;
    int face_info_size;
    fclaw2d_patch_relation_t prel;

    FCLAW_ASSERT (ind != NULL && ind->ready);
    FCLAW_ASSERT (domain == ind->domain);

    FCLAW_ASSERT (0 <= ghostno && ghostno < domain->num_ghost_patches);
    FCLAW_ASSERT (0 <= faceno && faceno < P4EST_FACES);

    /* check the type of neighbor situation */
    face_info_size = 2 + 2 * P4EST_HALF;
    pi = (int *) ind->e->ghost_data[ghostno] + face_info_size * faceno;
    indirect_match_face (pi, &grproc, &grblockno, &grpatchno, &grfaceno);
    *rblockno = *grblockno;
    FCLAW_ASSERT (0 <= *rblockno && *rblockno < domain->num_blocks);
    *rfaceno = *grfaceno & ~(3 << 26);
    FCLAW_ASSERT (0 <= *rfaceno && *rfaceno < P4EST_FACES * P4EST_HALF);
    if (!(*grfaceno & (3 << 26)))
    {
        if (grproc[0] == -1)
        {
            /* optimize for the most likely case */
            for (sf = 0; sf < P4EST_HALF; sf++)
            {
                FCLAW_ASSERT (grproc[sf] == -1 && grpatchno[sf] == -1);
                rproc[sf] = rpatchno[sf] = -1;
            }
            return FCLAW2D_PATCH_BOUNDARY;
        }
        else
        {
            prel = FCLAW2D_PATCH_SAMESIZE;
        }
    }
    else
    {
        if (*grfaceno & (1 << 26))
        {
            FCLAW_ASSERT (!(*grfaceno & (1 << 27)));
            prel = FCLAW2D_PATCH_HALFSIZE;
        }
        else
        {
            FCLAW_ASSERT (*grfaceno & (1 << 27));
            prel = FCLAW2D_PATCH_DOUBLESIZE;
        }
    }

    /* aslign the remaining output values */
    for (sf = 0; sf < P4EST_HALF; sf++)
    {
        rproc[sf] = grproc[sf];
        rpatchno[sf] = grpatchno[sf];
    }
#ifdef FCLAW_ENABLE_DEBUG
    if (rproc[0] != -1)
    {
        FCLAW_ASSERT (0 <= rproc[0] && rproc[0] < domain->mpisize);
        FCLAW_ASSERT (rproc[0] != domain->mpirank);
        FCLAW_ASSERT (0 <= rpatchno[0] &&
                      rpatchno[0] < domain->num_ghost_patches);
    }
    for (sf = 1; sf < P4EST_HALF; sf++)
    {
        if (prel == FCLAW2D_PATCH_HALFSIZE && rproc[sf] != -1)
        {
            FCLAW_ASSERT (0 <= rproc[sf] && rproc[sf] < domain->mpisize);
            FCLAW_ASSERT (rproc[sf] != domain->mpirank);
            FCLAW_ASSERT (0 <= rpatchno[sf] &&
                          rpatchno[sf] < domain->num_ghost_patches);
        }
    }
#endif

    /* and return */
    return prel;
}

fclaw2d_patch_relation_t
fclaw2d_domain_indirect_corner_neighbor (fclaw2d_domain_t * domain,
                                         fclaw2d_domain_indirect_t * ind,
                                         int ghostno, int cornerno,
                                         int *rproc, int *rblockno,
                                         int *rpatchno, int *rcornerno)
{
    int *pi;
    int *grproc, *grblockno, *grpatchno, *grcornerno;
    int face_info_size, corner_info_size;
    fclaw2d_patch_relation_t prel;

    FCLAW_ASSERT (ind != NULL && ind->ready);
    FCLAW_ASSERT (domain == ind->domain);

    FCLAW_ASSERT (0 <= ghostno && ghostno < domain->num_ghost_patches);
    FCLAW_ASSERT (0 <= cornerno && cornerno < P4EST_CHILDREN);

    /* check the type of neighbor situation */
    face_info_size = 2 + 2 * P4EST_HALF;
    corner_info_size = 4;
    pi = (int *) ind->e->ghost_data[ghostno] + face_info_size * P4EST_FACES
        + corner_info_size * cornerno;
    indirect_match_corner (pi, &grproc, &grblockno, &grpatchno, &grcornerno);
    *rblockno = *grblockno;
    FCLAW_ASSERT (0 <= *rblockno && *rblockno < domain->num_blocks);
    *rcornerno = *grcornerno & ~(3 << 26);
    FCLAW_ASSERT (0 <= *rcornerno && *rcornerno < P4EST_CHILDREN);
    if (!(*grcornerno & (3 << 26)))
    {
        if (grproc[0] == -1)
        {
            /* optimize for the most likely case */
            FCLAW_ASSERT (grproc[0] == -1 && grpatchno[0] == -1);
            *rproc = *rpatchno = -1;
            return FCLAW2D_PATCH_BOUNDARY;
        }
        else
        {
            prel = FCLAW2D_PATCH_SAMESIZE;
        }
    }
    else
    {
        if (*grcornerno & (1 << 26))
        {
            FCLAW_ASSERT (!(*grcornerno & (1 << 27)));
            prel = FCLAW2D_PATCH_HALFSIZE;
        }
        else
        {
            FCLAW_ASSERT (*grcornerno & (1 << 27));
            prel = FCLAW2D_PATCH_DOUBLESIZE;
        }
    }

    /* aslign the remaining output values */
    *rproc = grproc[0];
    *rpatchno = grpatchno[0];
#ifdef FCLAW_ENABLE_DEBUG
    if (*rproc != -1)
    {
        FCLAW_ASSERT (0 <= *rproc && *rproc < domain->mpisize);
        FCLAW_ASSERT (*rproc != domain->mpirank);
        FCLAW_ASSERT (0 <= *rpatchno &&
                      *rpatchno < domain->num_ghost_patches);
    }
#endif

    /* and return */
    return prel;
}


#ifdef P4_TO_P8
fclaw3d_patch_relation_t
fclaw3d_domain_indirect_edge_neighbors (fclaw3d_domain_t * domain,
                                        fclaw3d_domain_indirect_t * ind,
                                        int ghostno, int edgeno,
                                        int rproc[2], int *rblockno,
                                        int rpatchno[2],
                                        int *redgeno)
{
    int *pi;
    int *grproc, *grblockno, *grpatchno, *gredgeno;
    int se;
    int face_info_size, corner_info_size, edge_info_size;
    fclaw3d_patch_relation_t prel;

    FCLAW_ASSERT (ind != NULL && ind->ready);
    FCLAW_ASSERT (domain == ind->domain);

    FCLAW_ASSERT (0 <= ghostno && ghostno < domain->num_ghost_patches);
    FCLAW_ASSERT (0 <= edgeno && edgeno < P8EST_EDGES);

    /* check the type of neighbor situation */
    face_info_size = 2 + 2 * P8EST_HALF;
    corner_info_size = 4;
    edge_info_size = 2 + 2 * 2;
    pi = (int *) ind->e->ghost_data[ghostno] 
                    + face_info_size * P8EST_FACES
                    + corner_info_size * P8EST_CHILDREN
                    + edge_info_size * edgeno;
    indirect_match_edge (pi, &grproc, &grblockno, &grpatchno, &gredgeno);
    *rblockno = *grblockno;
    FCLAW_ASSERT (0 <= *rblockno && *rblockno < domain->num_blocks);
    *redgeno = *gredgeno & ~(3 << 26);
    FCLAW_ASSERT (0 <= *redgeno && *redgeno < P8EST_EDGES);
    if (!(*gredgeno & (3 << 26)))
    {
        if (grproc[0] == -1)
        {
            /* optimize for the most likely case */
            for (se = 0; se < 2; se++)
            {
                FCLAW_ASSERT (grproc[se] == -1 && grpatchno[se] == -1);
                rproc[se] = rpatchno[se] = -1;
            }
            return FCLAW3D_PATCH_BOUNDARY;
        }
        else
        {
            prel = FCLAW3D_PATCH_SAMESIZE;
        }
    }
    else
    {
        if (*gredgeno & (1 << 26))
        {
            FCLAW_ASSERT (!(*gredgeno & (1 << 27)));
            prel = FCLAW3D_PATCH_HALFSIZE;
        }
        else
        {
            FCLAW_ASSERT (*gredgeno & (1 << 27));
            prel = FCLAW3D_PATCH_DOUBLESIZE;
        }
    }

    /* aslign the remaining output values */
    for (se = 0; se < 2; se++)
    {
        rproc[se] = grproc[se];
        rpatchno[se] = grpatchno[se];
    }
#ifdef FCLAW_ENABLE_DEBUG
    if (rproc[0] != -1)
    {
        FCLAW_ASSERT (0 <= rproc[0] && rproc[0] < domain->mpisize);
        FCLAW_ASSERT (rproc[0] != domain->mpirank);
        FCLAW_ASSERT (0 <= rpatchno[0] &&
                      rpatchno[0] < domain->num_ghost_patches);
    }
    if (prel == FCLAW3D_PATCH_HALFSIZE && rproc[1] != -1)
    {
        FCLAW_ASSERT (0 <= rproc[1] && rproc[1] < domain->mpisize);
        FCLAW_ASSERT (rproc[1] != domain->mpirank);
        FCLAW_ASSERT (0 <= rpatchno[1] &&
                      rpatchno[1] < domain->num_ghost_patches);
    }
#endif

    /* and return */
    return prel;
}
#endif /* P4_TO_P8 */

void
fclaw2d_domain_indirect_destroy (fclaw2d_domain_t * domain,
                                 fclaw2d_domain_indirect_t * ind)
{
    FCLAW_ASSERT (ind != NULL && ind->ready);
    FCLAW_ASSERT (domain == ind->domain);

    fclaw2d_domain_free_after_exchange (domain, ind->e);
    FCLAW_FREE (ind);
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

int
fclaw2d_domain_is_meta (fclaw2d_domain_t * domain)
{
    FCLAW_ASSERT (domain != NULL);
    if (domain->local_num_patches == -1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void
fclaw2d_domain_init_meta (fclaw2d_domain_t * domain, int mpirank)
{
    FCLAW_ASSERT(domain != NULL);

    /* initialize to -1 and set pointers to NULL */
    memset (domain, -1, sizeof (fclaw2d_domain_t));
    domain->mpicomm = sc_MPI_COMM_NULL;
    domain->blocks = NULL;
    domain->exchange_patches = NULL;
    domain->ghost_patches = NULL;
    domain->mirror_target_levels = NULL;
    domain->ghost_target_levels = NULL;
    domain->pp = NULL;
    domain->attributes = NULL;
    domain->just_adapted = 0;
    domain->just_partitioned = 0;
    domain->pp_owned = 0;

    domain->local_num_patches = -1; /* mark as meta domain */

    domain->mpirank = mpirank; /* set mpirank to provide context information */
}

#if 0

/* This is how it used to work.  Never needed it */

/** Access a local or off-processor patch by its block and patch number.
 * \param [in] domain   Valid domain.
 * \param [in] blockno  For a local patch the number of the block.
 *                      To request a remote patch set this number to -1.
 * \param [in] patchno  For a local patch the number of the patch
 *                      relative to its block.  Otherwise the number
 *                      of the remote patch as returned by
 *                      \ref fclaw2d_patch_face_neighbors and friends.
 * \return              The patch that was requested.
 */
fclaw2d_patch_t *fclaw2d_domain_get_patch (fclaw2d_domain_t * domain,
                                           int blockno, int patchno);

#endif /* ! P4_TO_P8 */
