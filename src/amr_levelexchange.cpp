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
#include "amr_utils.H"

/* NOTE: Do we need the extern "C" here?  We're passing callbacks
   to C iterators but maybe C++ handles it just fine. */
#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

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

    // int numfaces = get_faces_per_patch(domain);

    for (int iface = 0; iface < NumFaces; iface++)
    {
        // Output arguments
        int neighbor_block_idx;
        int neighbor_patch_idx[p4est_refineFactor];
        int ref_flag;   // = -1, 0, 1
        int *ref_flag_ptr = &ref_flag;

        get_face_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iface,
                           &neighbor_block_idx,
                           neighbor_patch_idx,
                           &ref_flag_ptr);

        if (ref_flag_ptr == NULL)
        {
            // No neighbors - we are at a physical boundary
        }
        else if (ref_flag == 0)
        {
            // We have a neighbor patch at the same level
            fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
            fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[0]];
            ClawPatch *neighbor_cp = get_clawpatch(neighbor_patch);
            if (this_block_idx == neighbor_block_idx)
            {
                if (iface % 2 == 1)
                {
                    // Do high side exchange only
                    int idir = iface/2;   // this rounds down, right?  1/2 = 0; 3/2 = 1, etc.
                    // Exchange between 'this_patch' and 'neighbor patch(es)' in direction 'idir'
                    this_cp->exchange_face_ghost(idir,neighbor_cp);
                }
            }
            else
            {
                // Initiate exchange from block 0
                this_cp->mb_exchange_face_ghost(iface,neighbor_cp);
            }
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
            int corner_patch_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag;

            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                &corner_block_idx,
                                &corner_patch_idx,
                                &ref_flag_ptr);

            if (ref_flag_ptr == NULL)
            {
                // We should never get here, since we only call 'get_corner_neighbors' in a
                // situation in which we are sure we have neighbors.
            }
            else if (ref_flag == 0)
            {
                fclaw2d_block_t *corner_block = &domain->blocks[corner_block_idx];
                fclaw2d_patch_t *corner_patch = &corner_block->patches[corner_patch_idx];
                ClawPatch *corner_cp = get_clawpatch(corner_patch);
                if (this_block_idx == corner_block_idx)
                {
                    // Exchanging at the same level on the same block.
                    if (icorner % 2 == 1)
                    {
                        // Only initiate exchanges from high side corners when on the same block
                        this_cp->exchange_corner_ghost(icorner,corner_cp);
                    }
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

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

/* -------------------------------------------------------------------
   Main routine in this file
   ------------------------------------------------------------------- */
void level_exchange(fclaw2d_domain_t *domain, int level)
{
    fclaw2d_domain_iterate_level(domain, level,
                                 cb_level_face_exchange, (void *) NULL);

    // Do corner exchange only after physical boundary conditions have been set on all patches,
    // since corners may overlap phyical ghost cell region of neighboring patch.
    fclaw2d_domain_iterate_level(domain, level, cb_level_corner_exchange, (void *) NULL);
}
