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

#include "amr_forestclaw.H"


// ------------------------------------------------------------------------------
// Several routines that interact with p4est to get neighbor and corner information,
// and to determine if we are at physical boundaries or not.
// ------------------------------------------------------------------------------
void get_face_neighbors(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        int iside,
                        int *neighbor_block_idx,
                        int neighbor_patch_idx[],
                        int **ref_flag_ptr)
{
    const int p4est_refineFactor = get_p4est_refineFactor(domain);
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;

    fclaw2d_patch_relation_t neighbor_type =
        fclaw2d_patch_face_neighbors(domain,
                                     this_block_idx,
                                     this_patch_idx,
                                     iside,
                                     rproc,
                                     &rblockno,
                                     rpatchno,
                                     &rfaceno);


    // neighbor_type is one of :
    // FCLAW2D_PATCH_BOUNDARY,
    // FCLAW2D_PATCH_HALFSIZE,
    // FCLAW2D_PATCH_SAMESIZE,
    // FCLAW2D_PATCH_DOUBLESIZE

    *neighbor_block_idx = rblockno;

    if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
    {
        // Edge is a physical boundary
        // Set the pointer to NULL rather than come up with some bogus value for ref_flag.
        *ref_flag_ptr = NULL;
     }
    else
    {
        if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            // Neighbors are finer grids
            **ref_flag_ptr = 1;
            for(int ir = 0; ir < p4est_refineFactor; ir++)
            {
                neighbor_patch_idx[ir] = rpatchno[ir];
            }
        }
        else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            // Neighbor is at the same level
            **ref_flag_ptr = 0;
            neighbor_patch_idx[0] = rpatchno[0];
        }
        else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
        {
            // Neighbor is a coarser grid
            **ref_flag_ptr = -1;
            neighbor_patch_idx[0] = rpatchno[0];
        }
        else
        {
            // This didn't compile for me...
            // *(int *) 0 = 0;     // This must not happen
        }
    }
}

void get_corner_neighbor(fclaw2d_domain_t *domain,
                         int this_block_idx,
                         int this_patch_idx,
                         int icorner,
                         int *corner_block_idx,
                         int *corner_patch_idx,
                         int **ref_flag_ptr)
{
    fclaw2d_patch_relation_t neighbor_type;
    int rproc;
    int has_corner_neighbor =
      fclaw2d_patch_corner_neighbors (domain, this_block_idx, this_patch_idx,
        icorner, &rproc, corner_block_idx, corner_patch_idx, &neighbor_type);

    if (!has_corner_neighbor)
    {
        // printf("get_corner_neighbors : Patch %d at corner %d does not have any "
        //        "corner neighbors\n",this_patch_idx,icorner);
        *ref_flag_ptr = NULL;
    }
    else if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
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

void get_block_boundary(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        bool *intersects_block)
{
    const int p4est_refineFactor = get_p4est_refineFactor(domain);
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;
    const int numfaces = get_faces_per_patch(domain);

    for (int iside = 0; iside < numfaces; iside++)
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
        // neighbor_type is one of :
        // FCLAW2D_PATCH_BOUNDARY,
        // FCLAW2D_PATCH_HALFSIZE,
        // FCLAW2D_PATCH_SAMESIZE,
        // FCLAW2D_PATCH_DOUBLESIZE
        if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
        {
            // 'iside' is a physical boundary
            intersects_block[iside] = true;
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
                       bool *intersects_bc)
{
    const int numfaces = get_faces_per_patch(domain);
    int bdry[numfaces];
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);
    for(int i = 0; i < numfaces; i++)
    {
        // Physical boundary conditions
        intersects_bc[i] = bdry[i] == 1;
    }
}
