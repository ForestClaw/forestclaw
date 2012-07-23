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

const double fclaw2d_root_len = (double) P4EST_ROOT_LEN;
const double fclaw2d_smallest_h = 1. / (double) P4EST_ROOT_LEN;

void
fclaw2d_domain_iterate_level (fclaw2d_domain_t *domain, int level,
			      fclaw2d_patch_callback_t pcb, void *user)
{
  int			i, j;
  fclaw2d_block_t	*block;
  fclaw2d_patch_t	*patch;

  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  for (i = 0; i < domain->num_blocks; ++i) {
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

int
fclaw2d_patch_boundary_type (fclaw2d_domain_t *domain,
                int blockno, int patchno, int boundaries[P4EST_FACES])
{
  int                   faceno;
  int                   anyboundary;
  int8_t                qtf;
  p4est_locidx_t        totalleaf;
  p4est_locidx_t        qtq;
  p4est_t               *p4est;
  p4est_tree_t          *tree;
  p4est_mesh_t          *mesh;
  fclaw2d_block_t       *block;

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
  for (faceno = 0; faceno < P4EST_FACES; ++faceno) {
    qtq = mesh->quad_to_quad[P4EST_FACES * totalleaf + faceno];
    qtf = mesh->quad_to_face[P4EST_FACES * totalleaf + faceno];
    if (qtq == totalleaf && qtf == faceno) {
      anyboundary = 1;
      P4EST_ASSERT (block->mthbc[faceno] > 0);
      boundaries[faceno] = block->mthbc[faceno];
    }
    else {
      boundaries[faceno] = 0;
    }
  }

  return anyboundary;
}
