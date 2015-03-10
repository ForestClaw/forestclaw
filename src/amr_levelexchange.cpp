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

#include "fclaw2d_ghost_fill.h"
#include "amr_includes.H"

/* -----------------------------------------------------------------
   Exchange corner and face information at same level
   ----------------------------------------------------------------- */
#if 0
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
    transform_data.neighbor_patch = NULL;  // gets filled in below.


    // int numfaces = get_faces_per_patch(domain);

    for (int iface = 0; iface < NumFaces; iface++)
    {
        // Output arguments
        int neighbor_block_idx;
        int ref_flag;   // = -1, 0, 1
        int *ref_flag_ptr = &ref_flag;
        int fine_grid_pos;
        int *fine_grid_pos_ptr = &fine_grid_pos;

        /* Get the face neighbor relative to the neighbor's coordinate
           orientation (this isn't used here) */
        int iface_neighbor;
        int *iface_neighbor_ptr = &iface_neighbor;

        fclaw2d_patch_t* ghost_patches[p4est_refineFactor];


        transform_data.iface = iface;
        get_face_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iface,
                           &neighbor_block_idx,
                           ghost_patches,
                           &ref_flag_ptr,
                           &fine_grid_pos_ptr,
                           &iface_neighbor_ptr,
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

            transform_data.neighbor_patch = ghost_patches[0];
            this_cp->exchange_face_ghost(iface,neighbor_cp,&transform_data);

#if 0
            /* This is now done for all boundaries */
            if (this_block_idx == neighbor_block_idx)
            {
                transform_data.neighbor_patch = ghost_patches[0];
                this_cp->exchange_face_ghost(iface,neighbor_cp,&transform_data);
            }
            else
            {
                this_cp->mb_exchange_face_ghost(iface,neighbor_cp);
            }
#endif

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
    fclaw_bool intersects_bdry[NumFaces];
    fclaw_bool intersects_block[NumFaces];
    fclaw_bool is_block_corner;
    fclaw_bool is_interior_corner;

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bdry);

    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    /* Transform data needed at multi-block boundaries */
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.


    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        get_corner_type(domain,icorner,
                        intersects_bdry,
                        intersects_block,
                        &is_interior_corner,
                        &is_block_corner,
                        &transform_data.iface);


        if (is_interior_corner)  /* Interior to the domain, not necessarily to a block */
        {
            int corner_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag;
            int block_corner_count;
            fclaw2d_patch_t *ghost_patch;

            transform_data.icorner = icorner;
            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                transform_data.iface,
                                is_block_corner,
                                &corner_block_idx,
                                &ghost_patch,
                                &ref_flag_ptr,
                                &block_corner_count,
                                transform_data.transform);
            if (ref_flag_ptr == NULL)
            {
                /* Nothing happens here */
            }
            else if (ref_flag == 0)
            {
                ClawPatch *this_cp = get_clawpatch(this_patch);
                fclaw2d_patch_t* corner_patch = ghost_patch;
                ClawPatch *corner_cp = get_clawpatch(corner_patch);
                transform_data.neighbor_patch = corner_patch;
                if (!is_block_corner)
                {
                    this_cp->exchange_corner_ghost(icorner,corner_cp,
                                                   &transform_data);
                }
                else
                {
                    if (ispillowsphere_())
                    {
                        this_cp->mb_exchange_corner_ghost(icorner,intersects_block,
                                                          corner_cp,is_block_corner);
                    }
                }
            }
        }
    }
}
#endif

/* -------------------------------------------------------------------
   Main routine in this file
   ------------------------------------------------------------------- */
void level_exchange(fclaw2d_domain_t *domain, int level)
{
    fclaw2d_exchange_info_t filltype;
    filltype.exchange_type = FCLAW2D_COPY;
    filltype.grid_type = FCLAW2D_IS_COARSE;
    filltype.time_interp = fclaw_false;

    /* face exchanges */
    fclaw2d_domain_iterate_level(domain, level, cb_face_fill,
                                 (void *) &filltype);

    /* corner exchanges */
    fclaw2d_domain_iterate_level(domain, level, cb_corner_fill,
                                 (void *) &filltype);

}
