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

void cb_face_fill(fclaw2d_domain_t *domain,
                  fclaw2d_patch_t *this_patch,
                  int this_block_idx,
                  int this_patch_idx,
                  void *user)
{
    fclaw2d_exchange_info_t *filltype = (fclaw2d_exchange_info_t*) user;
    fclaw_bool time_interp = filltype->time_interp;
    fclaw_bool is_coarse = filltype->grid_type == FCLAW2D_IS_COARSE;
    fclaw_bool is_fine = filltype->grid_type == FCLAW2D_IS_FINE;

    fclaw_bool copy_from_neighbor = filltype->exchange_type == FCLAW2D_COPY;
    fclaw_bool average_from_neighbor = filltype->exchange_type == FCLAW2D_AVERAGE;
    fclaw_bool interpolate_to_neighbor = filltype->exchange_type == FCLAW2D_INTERPOLATE;


    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    /* Transform data needed at multi-block boundaries */
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.

    ClawPatch *this_cp = get_clawpatch(this_patch);
    for (int iface = 0; iface < NumFaces; iface++)
    {
        int idir = iface/2;

        /* Output arguments */
        int neighbor_block_idx;
        int relative_refinement_level;   // = -1, 0, 1
        int *ref_flag_ptr = &relative_refinement_level;
        int fine_grid_pos;
        int *fine_grid_pos_ptr = &fine_grid_pos;

        /* Get the face neighbor relative to the neighbor's coordinate
           orientation (this isn't used here) */
        int iface_neighbor;
        int *iface_neighbor_ptr = &iface_neighbor;

        fclaw2d_patch_t* neighbor_patches[p4est_refineFactor];

        transform_data.iface = iface;
        get_face_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iface,
                           &neighbor_block_idx,
                           neighbor_patches,
                           &ref_flag_ptr,
                           &fine_grid_pos_ptr,
                           &iface_neighbor_ptr,
                           transform_data.transform);

        fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
        if (ref_flag_ptr == NULL)
        {
            /* No face neighbors - we are at a physical boundary */
            continue;
        }

        if (copy_from_neighbor)
        {
            if (relative_refinement_level == 0)
            {
                /* We have a neighbor patch at the same level */
                fclaw2d_patch_t *neighbor_patch = neighbor_patches[0];
                ClawPatch *neighbor_cp = get_clawpatch(neighbor_patch);

                transform_data.neighbor_patch = neighbor_patches[0];
                this_cp->exchange_face_ghost(iface,neighbor_cp,&transform_data);
            }
        }
        else if (average_from_neighbor)
        {
            /* neighbors are at a finer level */
            if (relative_refinement_level == 1 && is_coarse)
            {
                /* Fill in ghost cells on 'this_patch' by averaging data from finer neighbors */
                fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
                for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                {
                    transform_data.neighbor_patch = neighbor_patches[igrid];
                    transform_data.fine_grid_pos = igrid;
                    ClawPatch* fine_neighbor_cp = get_clawpatch(neighbor_patches[igrid]);
                    this_cp->average_face_ghost(idir,iface,p4est_refineFactor,refratio,
                                                fine_neighbor_cp,time_interp,block_boundary,
                                                igrid,&transform_data);
                }
            }
            else if (relative_refinement_level == -1 && is_fine)
            {
                /* This happens only when running parallel */
                /* Found a fine grid with a parallel ghost patch as a neighbor */
                if (fclaw2d_patch_is_ghost (neighbor_patches[0]))
                {
                    ClawPatch *coarse_cp = get_clawpatch(neighbor_patches[0]);
                    ClawPatch *fine_cp = this_cp;
                    /* Figure out which grid we got */
                    int igrid = fine_grid_pos;

                    /* Swap out coarse and fine */
                    int iface_coarse = iface_neighbor;

                    fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
                    transform_data.this_patch = neighbor_patches[0];
                    transform_data.neighbor_patch = this_patch;
                    transform_data.fine_grid_pos = igrid;
                    coarse_cp->average_face_ghost(idir,iface_coarse,
                                                  p4est_refineFactor,refratio,
                                                  fine_cp,time_interp,block_boundary,
                                                  igrid, &transform_data);
                }
            }
        }
        else if (interpolate_to_neighbor)
        {
            if (relative_refinement_level == 1 && is_coarse)  // neighbors are at finer level
            {
                /* Fill in ghost cells on 'neighbor_patch' by interpolation */
                for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                {
                    ClawPatch *fine_neighbor_cp = get_clawpatch(neighbor_patches[igrid]);
                    transform_data.neighbor_patch = neighbor_patches[igrid];
                    transform_data.fine_grid_pos = igrid;
                    this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                    refratio,fine_neighbor_cp,
                                                    time_interp,
                                                    block_boundary,igrid,
                                                    &transform_data);
                }
            }
            else if (relative_refinement_level == -1 && is_fine)
            {
                /* Found a fine grid with a parallel coarse grid neighbor? */
                if (fclaw2d_patch_is_ghost (neighbor_patches[0]))
                {
                    ClawPatch *coarse_cp = get_clawpatch(neighbor_patches[0]);
                    ClawPatch *fine_cp = this_cp;
                    /* Figure out which grid we got */
                    int igrid = fine_grid_pos;

                    /* Swap out coarse and fine */
                    int iface_coarse = iface_neighbor;

                    /* This call can generate a floating point exception, for reasons I don't
                       completely understand, (er don't understand at all...) Valgrind doesn't
                       find any errors, and results do not have any NANs.  Strange
                    */
                    coarse_cp->interpolate_face_ghost(idir,iface_coarse,
                                                      p4est_refineFactor,refratio,
                                                      fine_cp,time_interp,
                                                      block_boundary,
                                                      igrid,
                                                      &transform_data);
                }
            }
        }
    }
}
