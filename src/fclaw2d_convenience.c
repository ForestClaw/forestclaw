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

static fclaw2d_domain_t		*
fclaw2d_domain_new (p4est_wrap_t *wrap, int mx, int my)
{
  int			i, j;
  int			face;
  int			nb;
  p4est_connectivity_t	*conn = wrap->conn;
  p4est_tree_t          *tree;
  p4est_quadrant_t      *quad;
  fclaw2d_domain_t	*domain;
  fclaw2d_block_t	*block;
  fclaw2d_patch_t	*patch;

  domain = P4EST_ALLOC_ZERO (fclaw2d_domain_t, 1);
  domain->mx = mx;
  domain->my = my;
  domain->pp = wrap;
  domain->num_blocks = nb = (int) conn->num_trees;
  domain->blocks = P4EST_ALLOC_ZERO (fclaw2d_block_t, domain->num_blocks);
  for (i = 0; i < nb; ++i) {
    block = domain->blocks + i;
    tree = p4est_tree_array_index (wrap->p4est->trees, (p4est_topidx_t) i);
    block->xlower = 0.;
    block->xupper = 1.;
    block->ylower = 0.;
    block->yupper = 1.;
    for (face = 0; face < P4EST_FACES; ++face) {
      if (conn->tree_to_tree[P4EST_FACES * nb + face] != (p4est_topidx_t) i ||
          conn->tree_to_face[P4EST_FACES * nb + face] != (int8_t) face) {
        block->mthbc[face] = 1;
      }
    }
    block->num_patches = (int) tree->quadrants.elem_count;
    block->patches = P4EST_ALLOC_ZERO (fclaw2d_patch_t, block->num_patches);
    for (j = 0; j < block->num_patches; ++j) {
      patch = block->patches + j;
      quad = p4est_quadrant_array_index (&tree->quadrants, (size_t) j);
      patch->level = (int) quad->level;
      patch->xlower = 0.;
      patch->xupper = 1.;
      patch->ylower = 0.;
      patch->yupper = 1.;
    }
  }

  return domain;
}

fclaw2d_domain_t		*
fclaw2d_domain_new_unitsquare (MPI_Comm mpicomm, int mx, int my)
{
  p4est_wrap_t		*wrap;

  wrap = p4est_wrap_new (mpicomm, 0);

  return fclaw2d_domain_new (wrap, mx, my);
}

void
fclaw2d_domain_destroy (fclaw2d_domain_t *domain)
{
  int			i;
  p4est_wrap_t		*wrap = domain->pp;
  fclaw2d_block_t	*block;

  for (i = 0; i < domain->num_blocks; ++i) {
    block = domain->blocks + i;
    P4EST_FREE (block->patches);
  }
  P4EST_FREE (domain->blocks);
  P4EST_FREE (domain);

  p4est_wrap_destroy (wrap);
}
