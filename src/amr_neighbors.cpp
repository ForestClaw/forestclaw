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
            /* 'iside' is a physical boundary (and so it must also be a
               block boundary */
            intersects_block[iside] = fclaw_true;
        }
        else
        {
            // We have a neighbor patch on block 'rblockno'.
            intersects_block[iside] = this_block_idx != rblockno;
        }
    }
}

// This is needed by other routines, so we don't set it to static.
void get_phys_boundary(fclaw2d_domain_t *domain,
                       int this_block_idx,
                       int this_patch_idx,
                       fclaw_bool *intersects_bdry)
{
    // const int numfaces = get_faces_per_patch(domain);
    int bdry[NumFaces];
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);
    for(int i = 0; i < NumFaces; i++)
    {
        // Physical boundary conditions
        intersects_bdry[i] = bdry[i] == 1;
    }
}
