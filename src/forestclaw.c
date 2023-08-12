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

#include <forestclaw.h>
#include <p4est.h>
#include <p8est.h>

void
fclaw_domain_attribute_add (fclaw_domain_t * domain,
                              const char *name, void *attribute)
{
    sc_keyvalue_t *a = domain->attributes;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (!sc_keyvalue_exists (a, name));
    sc_keyvalue_set_pointer (a, name, attribute);
}

void *
fclaw_domain_attribute_access (fclaw_domain_t * domain,
                                 const char *name, void *default_attr)
{
    sc_keyvalue_t *a = domain->attributes;

    FCLAW_ASSERT (a != NULL);
    return sc_keyvalue_get_pointer (a, name, default_attr);
}

void
fclaw_domain_attribute_remove (fclaw_domain_t * domain, const char *name)
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

void
fclaw_domain_iterate_level (fclaw_domain_t * domain, int level,
                              fclaw_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

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
fclaw_domain_iterate_patches (fclaw_domain_t * domain,
                                fclaw_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

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
fclaw_domain_iterate_families (fclaw_domain_t * domain,
                                 fclaw_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

    for (i = 0; i < domain->num_blocks; i++)
    {
        block = domain->blocks + i;
        for (j = 0; j < block->num_patches; j++)
        {
            patch = block->patches + j;
            if (fclaw_patch_is_first_sibling (patch))
            {
#ifdef FCLAW_ENABLE_DEBUG
                int k;
                for (k = 0; k < fclaw_domain_num_children(domain); ++k)
                {
                    FCLAW_ASSERT (j + k < block->num_patches);
                    FCLAW_ASSERT (fclaw_patch_childid (patch + k) == k);
                }
#endif
                pcb (domain, patch, i, j, user);
                j += fclaw_domain_num_children(domain) - 1;
            }
        }
    }
}

int
fclaw_domain_dimension (const fclaw_domain_t * domain)
{
    return domain->dim;
}

int
fclaw_domain_num_children (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_CHILDREN : P8EST_CHILDREN;
}

int
fclaw_domain_num_faces (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_FACES : P8EST_FACES;
}

int
fclaw3d_domain_num_edges (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? 0 : P8EST_EDGES;
}

int
fclaw_domain_num_corners (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_CHILDREN : P8EST_CHILDREN;
}

int
fclaw_domain_num_face_corners (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_HALF : P8EST_HALF;
}

int
fclaw_domain_num_orientations (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? (P4EST_FACES * P4EST_HALF) : (P8EST_FACES * P8EST_HALF);
}

void
fclaw_domain_corner_faces (const fclaw_domain_t * domain,
                             int icorner, int faces[P4EST_DIM])
{
    FCLAW_ASSERT (0 <= icorner && icorner < fclaw_domain_num_children(domain));
    if(domain->dim == 2)
    {
        faces[0] = p4est_corner_faces[icorner][0];
        faces[1] = p4est_corner_faces[icorner][1];
    }
    else
    {
        faces[0] = p8est_corner_faces[icorner][0];
        faces[1] = p8est_corner_faces[icorner][1];
        faces[2] = p8est_corner_faces[icorner][2];
    }
}

int
fclaw_patch_corner_dimension (const fclaw_patch_t * patch, int cornerno)
{
    const int childid = fclaw_patch_childid (patch);

    if(patch->dim == 2)
    {
        FCLAW_ASSERT (0 <= cornerno && cornerno < P4EST_CHILDREN);
        return (patch->level == 0 ||
                cornerno == childid ||
                cornerno == P4EST_CHILDREN - 1 - childid) ? 0 : 1;
    }
    else
    {
        FCLAW_ASSERT (0 <= cornerno && cornerno < P8EST_CHILDREN);
        return (patch->level == 0 ||
                cornerno == childid ||
                cornerno == P8EST_CHILDREN - 1 - childid) ? 0 : 1;
    }
}

int
fclaw_patch_childid (const fclaw_patch_t * patch)
{
    int childid;
    if(patch->dim == 2)
    {
        childid = patch->flags & FCLAW2D_PATCH_CHILDID;
        FCLAW_ASSERT (0 <= childid && childid < P4EST_CHILDREN);
    } 
    else 
    {
        childid = patch->flags & FCLAW3D_PATCH_CHILDID;
        FCLAW_ASSERT (0 <= childid && childid < P4EST_CHILDREN);
    }

    return childid;
}

int
fclaw_patch_is_first_sibling (const fclaw_patch_t * patch)
{
    if(patch->dim == 2)
    {
        return patch->flags & FCLAW2D_PATCH_FIRST_SIBLING ? 1 : 0;
    }
    else
    {
        return patch->flags & FCLAW3D_PATCH_FIRST_SIBLING ? 1 : 0;
    }
}

int
fclaw_patch_is_ghost (const fclaw_patch_t * patch)
{
    if(patch->dim == 2)
    {
        return patch->flags & FCLAW2D_PATCH_IS_GHOST ? 1 : 0;
    }
    else
    {
        return patch->flags & FCLAW3D_PATCH_IS_GHOST ? 1 : 0;
    }
}