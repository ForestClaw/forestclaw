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

#include "amr_includes.H"


// ------------------------------------------------------------------------------
// Several routines that interact with p4est to get neighbor and corner information,
// and to determine if we are at physical boundaries or not.
// ------------------------------------------------------------------------------

void get_face_neighbors(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        int iside,
                        int *neighbor_block_idx,
                        fclaw2d_patch_t* neighbor_patches[],
                        int **ref_flag_ptr,
                        int **fine_grid_pos_ptr,
                        int **iface_neighbor_ptr,
                        int ftransform[])
{
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;
    int num_neighbors;

    for(int ir = 0; ir < p4est_refineFactor; ir++)
    {
        neighbor_patches[ir] = NULL;
    }

    fclaw2d_patch_relation_t neighbor_type =
        fclaw2d_patch_face_neighbors(domain,
                                     this_block_idx,
                                     this_patch_idx,
                                     iside,
                                     rproc,
                                     &rblockno,
                                     rpatchno,
                                     &rfaceno);


    /* ------------------------------
      neighbor_type is one of :
      FCLAW2D_PATCH_BOUNDARY,
      FCLAW2D_PATCH_HALFSIZE,
      FCLAW2D_PATCH_SAMESIZE,
      FCLAW2D_PATCH_DOUBLESIZE
      ------------------------------- */

    *neighbor_block_idx = rblockno;

    if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
    {
        /* Edge is a physical boundary
           Set the pointer to NULL rather than come up with some bogus value
           for ref_flag and iface_neighbor */
        *ref_flag_ptr = NULL;
        *iface_neighbor_ptr = NULL;
    }
    else
    {
        // Get encoding of transforming a neighbor coordinate across a face
        fclaw2d_patch_face_transformation (iside, rfaceno, ftransform);
        if (this_block_idx == rblockno)
        {
            // If we are within one patch this is a special case
            ftransform[8] = 4;
        }

        if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            **ref_flag_ptr = 0;
            *fine_grid_pos_ptr = NULL;
            num_neighbors = 1;
        }
        else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
        {
            **ref_flag_ptr = -1;
            **fine_grid_pos_ptr = rproc[1];    // Special storage for fine grid info
            num_neighbors = 1;
        }
        else if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            /* Patch has two neighbors */
            **ref_flag_ptr = 1; /* patches are at one level finer */
            *fine_grid_pos_ptr = NULL;
            num_neighbors = p4est_refineFactor;
        }
        else
        {
            printf ("Illegal fclaw2d_patch_face_neighbors return value\n");
            exit (1);
        }

        for(int ir = 0; ir < num_neighbors; ir++)
        {
            fclaw2d_patch_t *neighbor;
            if (rproc[ir] == domain->mpirank)
            {
                /* neighbor patch is local */
                fclaw2d_block_t *neighbor_block = &domain->blocks[rblockno];
                neighbor = &neighbor_block->patches[rpatchno[ir]];
            }
            else
            {
                /* neighbor patch is on a remote processor */
                neighbor = &domain->ghost_patches[rpatchno[ir]];
            }
            neighbor_patches[ir] = neighbor;
        }
        **iface_neighbor_ptr = iside;
        fclaw2d_patch_face_swap(*iface_neighbor_ptr,&rfaceno);

    }
}

/* -----------------------------------------------
   Three corner cases to consider :
   (1) Corner is at a block corner. In this case,
       the number of blocks meeting at a corner
       is variable (2-pillow grid; 3-cubed sphere;
       4-brick; 5-???).  For now, this "is_block_corner"
       is handled as special cases.
   (2) Corner is at a block boundary, but not at a block
       corner.  Four patches meet at this corner, two of
       which will be on the remote block.  This corner
       has a face on 'this_patch' at iface, and at
       'rface' on the remote block.
   (3) Corner is internal to a single block.  This is the
       easiest case to handle.
   (4) Corner is a hanging node.  Don't need to do anything.

   In cases 2 and 3, a 'transform' between corner-adjacent
   patches is well defined.  In case 1, only when four
   patches meet is the transform well defined.

   Situations in which we might not have a corner
   adjacent neighbor include : cubed sphere case or
   hanging node case.
   ----------------------------------------------- */


void get_corner_neighbor(fclaw2d_domain_t *domain,
                         int this_block_idx,
                         int this_patch_idx,
                         int icorner,
                         int iface,
                         int *corner_block_idx,
                         fclaw2d_patch_t** ghost_patch,
                         int **ref_flag_ptr,
                         fclaw_bool is_block_corner,
                         int ftransform[])
{
    int rproc_corner;
    int corner_patch_idx;
    fclaw2d_patch_relation_t neighbor_type;
    fclaw_bool has_corner_neighbor = fclaw_false;

    if (is_block_corner)
    {
        /* Case (1) : 2, 3, or more  blocks meeting at a block corner */
        if (ispillowsphere_())
        {
            /* 2 patches meet at a corner.  This is handled as a special
               case */
            has_corner_neighbor = fclaw_true;
            int rpatchno[p4est_refineFactor];
            int rproc[p4est_refineFactor];
            int rfaceno;

            /* Use only faces 0 or 1 to get block data. */
            int iface = icorner % 2;
            neighbor_type =
                fclaw2d_patch_face_neighbors(domain,
                                             this_block_idx,
                                             this_patch_idx,
                                             iface,
                                             rproc,
                                             corner_block_idx,
                                             rpatchno,
                                             &rfaceno);

            int igrid;
            if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
            {
                /* igrid = 0 at corners 0,1 and (R-1) at corners 2,3,
                   where R = refinement factor */
                igrid = (icorner/2)*(p4est_refineFactor - 1);
            }
            else
            {
                /* Same size or double size grids. */
                igrid = 0;
            }

            corner_patch_idx = rpatchno[igrid];
            rproc_corner = rproc[igrid];
        }
        else
        {
            int rblockno;
            has_corner_neighbor =
                fclaw2d_patch_corner_neighbors(domain, this_block_idx, this_patch_idx,
                                               icorner, &rproc_corner, &rblockno,
                                               &corner_patch_idx, &neighbor_type);
            if (!has_corner_neighbor)
            {
                /* Case (1) : 3 patches meet at a corner, e.g. the
                   cubed sphere.   Nothing to do here */
            }
            else
            {
                /* Case (1) : 4 or more patches meet at a corner */
                printf("get_corner_neighbors (amr_neighbors.f) : " \
                       "Four or more patches meet at a block corner; " \
                       "(not implemented yet).\n");
                exit(0);
            }
        }
    }
    else
    {
        /* Cases (2) and (3) : We are either at a block edge or at an
           internal corner.  Unless we are a hanging node, we return
           a "transform" vector here. */
        int rblockno;
        has_corner_neighbor =
            fclaw2d_patch_corner_neighbors(domain, this_block_idx, this_patch_idx,
                                           icorner, &rproc_corner, &rblockno,
                                           &corner_patch_idx, &neighbor_type);

        if (!has_corner_neighbor)
        {
            /* 'icorner' is a hanging node.  See see below. */
        }
        else
        {
            *corner_block_idx = rblockno;
            if (this_block_idx != rblockno)
            {
                /* Case (2).  The corner is on a block edge (but is not a block corner) */

                /* Need to get remote face of remote corner grid.  We will
                   ignore most other returns from this function, including
                   rproc[], rblockno (== corner_block_idx, from above),
                   rpatchno and neighbor_return_type */
                int rfaceno;
                int rproc[p4est_refineFactor];
                int rpatchno;
                fclaw2d_patch_face_neighbors(domain,
                                             this_block_idx,
                                             this_patch_idx,
                                             iface,
                                             rproc,
                                             &rblockno,
                                             &rpatchno,
                                             &rfaceno);

                /* Get encoding of transforming a neighbor coordinate across a face */
                fclaw2d_patch_face_transformation (iface, rfaceno, ftransform);
            }
            else
            {
                /* Case (3) : Both patches are in the same block */
                ftransform[8] = 4;
            }
        }
    }

    /* ---------------------------------------------------------------------
       We now have a neighbor and possibly a transform.  Just need to figure
       out what level it is at.
       --------------------------------------------------------------------- */

    if (!has_corner_neighbor)
    {
        /* This can happen in the hanging node case, or in the cubed sphere case, at
           block corners */
        *ref_flag_ptr = NULL;
        *ghost_patch = NULL;
    }
    else
    {
        if (domain->mpirank != rproc_corner)
        {
            *ghost_patch = &domain->ghost_patches[corner_patch_idx];
        }
        else
        {
            fclaw2d_block_t *neighbor_block = &domain->blocks[*corner_block_idx];
            *ghost_patch = &neighbor_block->patches[corner_patch_idx];
        }

        /* ---------------------------
           neighbor_type is one of :
           FCLAW2D_PATCH_BOUNDARY,
           FCLAW2D_PATCH_HALFSIZE,
           FCLAW2D_PATCH_SAMESIZE,
           FCLAW2D_PATCH_DOUBLESIZE
           --------------------------- */

        if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            **ref_flag_ptr = 1;
        }
        else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            **ref_flag_ptr = 0;
        }
        else
        {
            **ref_flag_ptr = -1;
        }
    }
}

void get_block_boundary(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        fclaw_bool *intersects_block)
{
    // const int p4est_refineFactor = get_p4est_refineFactor(domain);
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;
    int ftransform[9];
    // const int numfaces = get_faces_per_patch(domain);

    for (int iside = 0; iside < NumFaces; iside++)
    {
        fclaw2d_patch_relation_t neighbor_type =
            fclaw2d_patch_face_neighbors(domain,
                                         this_block_idx,
                                         this_patch_idx,
                                         iside,
                                         rproc,
                                         &rblockno,
                                         rpatchno,
                                         &rfaceno);
        /* --------------------------
           neighbor_type is one of :
           FCLAW2D_PATCH_BOUNDARY,
           FCLAW2D_PATCH_HALFSIZE,
           FCLAW2D_PATCH_SAMESIZE,
           FCLAW2D_PATCH_DOUBLESIZE
           --------------------------*/

        if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
        {
            /* 'iside' is a physical boundary */
            intersects_block[iside] = true;
        }
        else
        {
            // This is demonstration code to infer the patch transformation.
            // We can do this at a later point when we really need to transform
            // a patch into it's face neighbor's coordinate system.
            // This is nontrivial only across a block boundary
            if (this_block_idx != rblockno) {
                fclaw2d_patch_face_transformation (iside, rfaceno, ftransform);
            }

            // We have a neighbor patch on block 'rblockno'.
            intersects_block[iside] = this_block_idx != rblockno;
        }
    }
}

// This is needed by other routines, so we don't set it to static.
void get_phys_boundary(fclaw2d_domain_t *domain,
                       int this_block_idx,
                       int this_patch_idx,
                       fclaw_bool *intersects_bc)
{
    // const int numfaces = get_faces_per_patch(domain);
    int bdry[NumFaces];
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);
    for(int i = 0; i < NumFaces; i++)
    {
        // Physical boundary conditions
        intersects_bc[i] = bdry[i] == 1;
    }
}
