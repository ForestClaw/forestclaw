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

/* -----------------------------------------------------------------
   Exchange corner and face information at same level
   ----------------------------------------------------------------- */

static
void cb_level_face_exchange(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            void *user)
{
    // const int p4est_refineFactor = get_p4est_refineFactor(domain);
    ClawPatch *this_cp = get_clawpatch(this_patch);

    // Transform data needed at multi-block boundaries
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    // The rest gets filled in below.


    // int numfaces = get_faces_per_patch(domain);

    for (int iface = 0; iface < NumFaces; iface++)
    {
        // Output arguments
        int neighbor_block_idx;
        int ref_flag;   // = -1, 0, 1
        int *ref_flag_ptr = &ref_flag;
        int fine_grid_pos;
        int *fine_grid_pos_ptr = &fine_grid_pos;
        fclaw2d_patch_t* ghost_patches[p4est_refineFactor];
        // int ftransform[9];


        transform_data.iface = iface;
        get_face_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iface,
                           &neighbor_block_idx,
                           ghost_patches,
                           &ref_flag_ptr,
                           &fine_grid_pos_ptr,
                           transform_data.transform);

        if (ref_flag_ptr == NULL)
        {
            // No neighbors - we are at a physical boundary
        }
        else if (ref_flag == 0)
        {
            /* We have a neighbor patch at the same level */
            fclaw2d_patch_t *neighbor_patch = ghost_patches[0];
            ClawPatch *neighbor_cp = get_clawpatch(neighbor_patch);

            /* This is now done for all boundaries */
            transform_data.neighbor_patch = ghost_patches[0];
            fclaw_cptr cptr;
            if (sizeof(cptr) != sizeof(&transform_data))
            {
                printf("amr_level_exchange : assumed size of ptr is incorrect; \
                        sizeof(cptr) = %u but sizeof(&transform_data) = %u\n",
                       sizeof(cptr), sizeof(&transform_data));
                exit(0);
            }
            this_cp->exchange_face_ghost(iface,neighbor_cp,(fclaw_cptr) &transform_data);

            /*
            if (this_block_idx == neighbor_block_idx)
            {
                // Exchange between 'this_patch' and 'neighbor patch(es)' at face 'iface'
                this_cp->exchange_face_ghost(iface,neighbor_cp);
            }
            else
            {
                // Initiate exchange from block 0
                this_cp->mb_exchange_face_ghost(iface,neighbor_cp);
            }
            */
        } /* Check return from neighbor */
    } /* loop over all faces */
}

static
void cb_level_corner_exchange(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              void *user)
{

    // const int numfaces = get_faces_per_patch(domain);
    fclaw_bool intersects_bc[NumFaces];
    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    fclaw_bool intersects_block[NumFaces];
    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    // Number of patch corners, not the number of corners in the domain!
    // const int numcorners = get_corners_per_patch(domain);


    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        // p4est has tons of lookup table like this, can be exposed similarly
        int corner_faces[SpaceDim];
        fclaw2d_domain_corner_faces(domain, icorner, corner_faces);

        // Both faces are at a physical boundary
        fclaw_bool is_phys_corner =
                intersects_bc[corner_faces[0]] && intersects_bc[corner_faces[1]];

        // Corner lies in interior of physical boundary edge.
        fclaw_bool corner_on_phys_face = !is_phys_corner &&
                (intersects_bc[corner_faces[0]] || intersects_bc[corner_faces[1]]);

        // Either a corner at a block boundary (but not a physical boundary),
        // or internal to a block
        fclaw_bool interior_corner = !corner_on_phys_face && !is_phys_corner;

        ClawPatch *this_cp = get_clawpatch(this_patch);
        if (is_phys_corner)
        {
            // We don't have to worry about this now.  It is
            // now taken care of by smart sweeping in the bc2 routine.
        }
        else if (corner_on_phys_face)
        {
            // Again, smart sweeping on in the bc2 routine should take care of these
            // corner cells.
        }
        else if (interior_corner)
        {
            // Both faces are at a block boundary
            fclaw_bool is_block_corner =
                intersects_block[corner_faces[0]] && intersects_block[corner_faces[1]];

            // We know corner 'icorner' has an adjacent patch.
            int corner_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag;
            fclaw2d_patch_t *ghost_patch;

            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                &corner_block_idx,
                                &ghost_patch,
                                &ref_flag_ptr,
                                is_block_corner);

            if (ref_flag_ptr == NULL)
            {
                // We should never get here, since we only call 'get_corner_neighbors' in a
                // situation in which we are sure we have neighbors.
            }
            else if (ref_flag == 0)
            {
                fclaw2d_patch_t* corner_patch = ghost_patch;
                ClawPatch *corner_cp = get_clawpatch(corner_patch);
                if (this_block_idx == corner_block_idx)
                {
                    this_cp->exchange_corner_ghost(icorner,corner_cp);
                }
                else
                {
                    // We are doing a corner exchange across blocks
                    this_cp->mb_exchange_corner_ghost(icorner,intersects_block,
                                                      corner_cp,is_block_corner);
                }
            }
        }
    }
}


/* -------------------------------------------------------------------
   Main routine in this file
   ------------------------------------------------------------------- */
void level_exchange(fclaw2d_domain_t *domain, int level)
{
    // Start exchanging

    /* face exchanges */
    fclaw2d_domain_iterate_level(domain, level,
                                 cb_level_face_exchange, (void *) NULL);

    /* Do corner exchange only after physical boundary conditions have
       been set on all patches, since corners may overlap physical ghost
       cell region of neighboring patch. ??? (where am I doing set_physbc?)
    */
    fclaw2d_domain_iterate_level(domain, level, cb_level_corner_exchange,
                                 (void *) NULL);

}
