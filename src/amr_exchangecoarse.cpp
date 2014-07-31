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

/* ******************************************************************************
   This file contains all the routines (averaging and interpolation between faces
   and corners) needed for a level to do an exchange (averaging and interpolation)
   with a coarser level.  The two main routines in this file are :

          void exchange_with_coarse(fclaw2d_domain_t *domain, const int& a_level);

   and the analogous routine, in the case that ghost cells are to be obtained from a
   temporary time interpolated grid :

       void exchange_with_coarse_time_interp(fclaw2d_domain_t *domain, const int& level,
                                             const int& coarse_step, const int& fine_step,
                                             const int& refratio)

   where "level" is the finer grid which should exchange with "level-1" (either
   at the same time at 'a_level', or at a time interpolated value.

   Routines in this file are :

      cb_corner_exchange, cb_corner_interpolate
      cb_face_average, cb_face_interpolate

   This is called from "advance_level" and it is assumed that "a_level" is not
   the finest grid.
   ******************************************************************************
*/

struct exchange_info
{
    fclaw_bool time_interp;
    fclaw_bool is_coarse;
    fclaw_bool is_fine;
    int level;
};


static
void cb_corner_average(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    exchange_info *e_info = (exchange_info*) user;
    fclaw_bool time_interp = e_info->time_interp;
    fclaw_bool is_coarse = e_info->is_coarse;
    fclaw_bool is_fine = e_info->is_fine;


    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;
    fclaw_bool intersects_bc[NumFaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    fclaw_bool intersects_block[NumFaces];
    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    // Number of patch corners, not the number of corners in the domain!
    // const int numcorners = get_corners_per_patch(domain);

    // Transform data needed at multi-block boundaries
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.

    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        // Get faces that intersect 'icorner'
        // There must be a clever way to do this...
        // p4est has tons of lookup table like this, can be exposed similarly
        int faces[SpaceDim];
        fclaw2d_domain_corner_faces(domain, icorner, faces);

        // Both faces are at a physical boundary
        fclaw_bool is_phys_corner =
                intersects_bc[faces[0]] && intersects_bc[faces[1]];

        // Corner lies in interior of physical boundary edge.
        fclaw_bool corner_on_phys_face = !is_phys_corner &&
                (intersects_bc[faces[0]] || intersects_bc[faces[1]]);

        // This corner may still be on a block boundary.
        fclaw_bool interior_corner = !corner_on_phys_face && !is_phys_corner;

        ClawPatch *this_cp = get_clawpatch(this_patch);
        if (is_phys_corner)
        {
            // We don't have to worry about this now.  It is
            // now taken care of by applying physical boundary conditions,
            // first in x, then in y.
        }
        else if (corner_on_phys_face)
        {
            // Taken care of by applying boundary conditions _after_ doing
            // face exchanges or coarse/fine interpolation at a face that
            // intersects (perpendicularly) the physical boundary.
        }
        else if (interior_corner)
        {
            // Both faces are at a block boundary
            fclaw_bool is_block_corner =
                intersects_block[faces[0]] && intersects_block[faces[1]];

            int corner_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag;
            fclaw2d_patch_t *neighbor_patch;

            transform_data.icorner = icorner;
            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                &corner_block_idx,
                                &neighbor_patch,
                                &ref_flag_ptr,
                                is_block_corner,
                                transform_data.transform);

            /* -------------------------------------
               Possible returns from ref_flag :
               FCLAW2D_PATCH_BOUNDARY : ref_flag_ptr = NULL
               FCLAW2D_PATCH_HALFSIZE : ref_flag = 1 (one level finer)
               FCLAW2D_PATCH_SAMESIZE : ref_flag = 0 (same level)
               FCLAW2D_PATCH_DOUBLESIZE : ref_flag = -1 (one level coarser)
               ------------------------------------- */

            if (ref_flag_ptr == NULL)
            {
                // No corner neighbor
            }
            else if (ref_flag == 1 && is_coarse)
            {
                /* Corner neighbor at a finer level, and so we need to average
                   that corner onto the coarser corner ghost cells */
                ClawPatch *corner_cp = get_clawpatch(neighbor_patch);
                transform_data.neighbor_patch = neighbor_patch;

                if (this_block_idx == corner_block_idx)
                {
                    this_cp->average_corner_ghost(icorner,refratio,corner_cp,
                                                  time_interp,
                                                  (fclaw_cptr) &transform_data);
                }
                else
                {
                    this_cp->mb_average_corner_ghost(icorner,refratio,corner_cp,time_interp,
                                                     is_block_corner,intersects_block);
                }
            }
            else if (ref_flag == -1 && is_fine)
            {
                /* Corner is at at coarser level.  If the corner patch is a parallel
                   ghost patch, then we should average onto it as well (even if it
                   is technically read only)  Otherwise, it won't have good data for
                   doing interpolation to this fine grid later */
                if (fclaw2d_patch_is_ghost (neighbor_patch))
                {
                    /* Corner is a parallel ghost patch.  Do the same as above, but
                       but from the perspective of the ghost (coarser) patch.*/
                    ClawPatch *corner_cp = this_cp;
                    ClawPatch *coarse_cp = get_clawpatch(neighbor_patch);

                    if (this_block_idx == corner_block_idx)
                    {
                        int icorner_coarse = 3-icorner;
                        coarse_cp->average_corner_ghost(icorner_coarse,
                                                        refratio,corner_cp,
                                                        time_interp,
                                                        (fclaw_cptr) &transform_data);
                    }
                    else
                    {
                        int icorner_coarse;
                        if (is_block_corner)
                        {
                            icorner_coarse = icorner;
                        }
                        else
                        {
                            if (icorner == 0)
                            {
                                if (intersects_block[0])
                                    icorner_coarse = 2;
                                else if (intersects_block[2])
                                    icorner_coarse = 1;
                            }
                            else if (icorner == 1)
                            {
                                if (intersects_block[1])
                                    icorner_coarse = 3;
                                else if (intersects_block[2])
                                    icorner_coarse = 0;
                            }
                            else if (icorner == 2)
                            {
                                if (intersects_block[0])
                                    icorner_coarse = 0;
                                else if (intersects_block[3])
                                    icorner_coarse = 3;
                            }
                            else if (icorner == 3)
                            {
                                if (intersects_block[1])
                                    icorner_coarse = 1;
                                else if (intersects_block[3])
                                    icorner_coarse = 2;
                            }
                        }
                        coarse_cp->mb_average_corner_ghost(icorner_coarse,
                                                           refratio,
                                                           corner_cp,
                                                           time_interp,
                                                           is_block_corner,
                                                           intersects_block);
                    }
                }
            }
        }
    }
}

static
void cb_corner_interpolate(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           void *user)
{
    exchange_info *e_info = (exchange_info*) user;
    fclaw_bool time_interp = e_info->time_interp;
    fclaw_bool is_fine = e_info->is_fine;

    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;
    fclaw_bool intersects_bc[NumFaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    fclaw_bool intersects_block[NumFaces];
    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    // Number of patch corners, not the number of corners in the domain!
    // const int numcorners = get_corners_per_patch(domain);

    // Transform data needed at multi-block boundaries
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.


    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        // Get faces that intersect 'icorner'
        // There must be a clever way to do this...
        // p4est has tons of lookup table like this, can be exposed similarly
        int faces[SpaceDim];
        fclaw2d_domain_corner_faces(domain, icorner, faces);

        // Both faces are at a physical boundary
        fclaw_bool is_phys_corner =
                intersects_bc[faces[0]] && intersects_bc[faces[1]];

        // Corner lies in interior of physical boundary edge.
        fclaw_bool corner_on_phys_face = !is_phys_corner &&
                (intersects_bc[faces[0]] || intersects_bc[faces[1]]);

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
                intersects_block[faces[0]] && intersects_block[faces[1]];

            int corner_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag;
            fclaw2d_patch_t *neighbor_patch;

            transform_data.icorner = icorner;
            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                &corner_block_idx,
                                &neighbor_patch,
                                &ref_flag_ptr,
                                is_block_corner,
                                transform_data.transform);

            if (ref_flag_ptr == NULL)
            {
                // no corner neighbor
            }
            else if (ref_flag == 1)
            {
                ClawPatch *corner_cp = get_clawpatch(neighbor_patch);

                transform_data.neighbor_patch = neighbor_patch;
                if (this_block_idx == corner_block_idx)
                {
                    this_cp->interpolate_corner_ghost(icorner,refratio,corner_cp,
                                                      time_interp,
                                                      (fclaw_cptr) &transform_data);
                }
                else
                {
                    this_cp->mb_interpolate_corner_ghost(icorner,refratio,
                                                         corner_cp,time_interp,
                                                         is_block_corner, intersects_block);
                }
            }
            else if (ref_flag == -1 && is_fine)
            {
	        /* Corner is at a coarser level.  Use this corner patch to do the interpolation */
                if (fclaw2d_patch_is_ghost(neighbor_patch))
                {
                    /* Corner is a parallel ghost patch.  Do the same as above, but
                       but from the perspective of the ghost (coarser) patch.*/
                    ClawPatch *corner_cp = this_cp;
                    ClawPatch *coarse_cp = get_clawpatch(neighbor_patch);

                    set_debug_info_(this_block_idx, this_patch_idx,this_patch->level);

                    if (this_block_idx == corner_block_idx)
                    {
                        int icorner_coarse = 3-icorner;
                        transform_data.icorner = icorner_coarse;
                        transform_data.this_patch = neighbor_patch;
                        transform_data.neighbor_patch = this_patch;
                        coarse_cp->interpolate_corner_ghost(icorner_coarse,
                                                            refratio,corner_cp,
                                                            time_interp,
                                                            (fclaw_cptr) &transform_data);
                    }
                    else
                    {
                        /* Works only for the sphere grid */
                        int icorner_coarse;
                        if (is_block_corner)
                        {
                            icorner_coarse = icorner;
                        }
                        else
                        {
                            if (icorner == 0)
                            {
                                if (intersects_block[0])
                                    icorner_coarse = 2;
                                else if (intersects_block[2])
                                    icorner_coarse = 1;
                            }
                            else if (icorner == 1)
                            {
                                if (intersects_block[1])
                                    icorner_coarse = 3;
                                else if (intersects_block[2])
                                    icorner_coarse = 0;
                            }
                            else if (icorner == 2)
                            {
                                if (intersects_block[0])
                                    icorner_coarse = 0;
                                else if (intersects_block[3])
                                    icorner_coarse = 3;
                            }
                            else if (icorner == 3)
                            {
                                if (intersects_block[1])
                                    icorner_coarse = 1;
                                else if (intersects_block[3])
                                    icorner_coarse = 2;
                            }
                        }
                        coarse_cp->mb_interpolate_corner_ghost(icorner_coarse,
                                                               refratio,corner_cp,time_interp,
                                                               is_block_corner,intersects_block);
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------
// Multi-level ghost cell operations
//   -- average fine grid ghost cells to coarse grid
//   -- interpolate coarse grid to fine grid ghost cells
// -----------------------------------------------------------------
static
void cb_face_average(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    exchange_info *e_info = (exchange_info*) user;
    fclaw_bool time_interp = e_info->time_interp;

    /* is_coarse == true : We are iterating over the coarser level -
                    look for fine grid neighbors
       is_fine == true : We are iterating over fine grids - look for
                         coarse grid off-proc ghost patches that we may need
                         to average to.
    */
    fclaw_bool is_coarse = e_info->is_coarse;
    fclaw_bool is_fine = e_info->is_fine;


    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    // Transform data needed at multi-block boundaries
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.

    ClawPatch *this_cp = get_clawpatch(this_patch);
    for (int idir = 0; idir < 2; idir++)
    {
        // Loop over low side and high side
        for (int iface = 2*idir; iface <= 2*idir + 1; iface++)
        {
            int neighbor_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag; // = -1, 0, 1
            int fine_grid_pos;
            int *fine_grid_pos_ptr = &fine_grid_pos;

            /* Get the face neighbor relative to the neighbor's coordinate
               orientation */
            int iface_neighbor;
            int *iface_neighbor_ptr = &iface_neighbor;

            fclaw2d_patch_t *neighbor_patches[p4est_refineFactor];

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

            if (ref_flag_ptr == NULL)
            {
                // no face neighbor
            }
            else if (ref_flag == 1 && is_coarse)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'this_patch' by averaging data from finer neighbors
                fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
                for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                {
                    transform_data.neighbor_patch = neighbor_patches[igrid];
                    transform_data.fine_grid_pos = igrid;
                    ClawPatch* fine_neighbor_cp = get_clawpatch(neighbor_patches[igrid]);
                    this_cp->average_face_ghost(idir,iface,p4est_refineFactor,refratio,
                                                fine_neighbor_cp,time_interp,block_boundary,
                                                igrid,(fclaw_cptr) &transform_data);
                }
            }
            else if (ref_flag == -1 && is_fine)
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

                    fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
#if 0
                    int iface_coarse;
                    if (!block_boundary)
                    {
                        if (iface == 0)
                            iface_coarse = 1;
                        else if (iface == 1)
                            iface_coarse = 0;
                        else if (iface == 2)
                            iface_coarse = 3;
                        else if (iface == 3)
                            iface_coarse = 2;
                    }
                    else
                    {
                        /* This only works for the sphere grid */
                        iface_coarse = iface;
                    }
#endif


                    transform_data.this_patch = neighbor_patches[0];
                    transform_data.neighbor_patch = this_patch;
                    transform_data.fine_grid_pos = igrid;
                    coarse_cp->average_face_ghost(idir,iface_coarse,p4est_refineFactor,refratio,
                                                  fine_cp,time_interp,block_boundary,
                                                  igrid, (fclaw_cptr) &transform_data);

                }

            } /* loop sides (iside = 0,1,2,3) */
        } /* loop over directions (idir = 0,1) */
    }
}



// Iterator over patches looking for finer neighbors
static
void cb_face_interpolate(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    exchange_info *e_info = (exchange_info*) user;
    fclaw_bool time_interp = e_info->time_interp;
    fclaw_bool is_coarse = e_info->is_coarse;
    fclaw_bool is_fine = e_info->is_fine;

    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    // Transform data needed at multi-block boundaries
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.


    ClawPatch *this_cp = get_clawpatch(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iface = 2*idir; iface <= 2*idir + 1; iface++)
        {
            int neighbor_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag; // = -1, 0, 1
            int fine_grid_pos;
            int *fine_grid_pos_ptr = &fine_grid_pos;

            /* Get the face neighbor relative to the neighbor's coordinate
               orientation */
            int iface_neighbor;
            int *iface_neighbor_ptr = &iface_neighbor;

            fclaw2d_patch_t *neighbor_patches[p4est_refineFactor];

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
                // no face neighbor - physical boundary?
            }
            else if (ref_flag == 1 && is_coarse)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'neighbor_patch' by interpolation
                for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                {
                    ClawPatch *fine_neighbor_cp = get_clawpatch(neighbor_patches[igrid]);
                    transform_data.neighbor_patch = neighbor_patches[igrid];
                    transform_data.fine_grid_pos = igrid;
                    this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                    refratio,fine_neighbor_cp,time_interp,
                                                    block_boundary,igrid,
                                                    (fclaw_cptr) &transform_data);
                }
            }
            else if (ref_flag == -1 && is_fine)
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
#if 0
                    int iface_coarse;
                    if (!block_boundary)
                    {
                        if (iface == 0)
                            iface_coarse = 1;
                        else if (iface == 1)
                            iface_coarse = 0;
                        else if (iface == 2)
                            iface_coarse = 3;
                        else if (iface == 3)
                            iface_coarse = 2;
                    }
                    else
                    {
                        /* This only works for the sphere grid */
                        iface_coarse = iface;
                    }
#endif

                    /* This call can generate a floating point exception, for reasons I don't
                       completely understand, (er don't understand at all...) Valgrind doesn't
                       find any errors, and results do not have any NANs.  Strange
                    */
                    coarse_cp->interpolate_face_ghost(idir,iface_coarse,
                                                      p4est_refineFactor,refratio,
                                                      fine_cp,time_interp,
                                                      block_boundary,
                                                      igrid,
                                                      (fclaw_cptr) &transform_data);
                }
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}


#if 0
static
void cb_setup_time_interp(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    /* This is called for all patches on the coarse level */
    ClawPatch *cp = get_clawpatch(this_patch);
    double &alpha = *((double*) user);

    /* This constructs a time interpolated version of the data on
       the coarser grid */
    cp->setup_for_time_interpolation(alpha);
}
#endif



/* ----------------------------------------------------------------------
   Main routine in this file.  This file assumes that both coarse and
   fine grids have valid interior data;  they only need to exchange (
   via interpolating and averaging) ghost cell values.
   -------------------------------------------------------------------- */

void exchange_with_coarse(fclaw2d_domain_t *domain,
                          int level, double t_level,
                          fclaw_bool time_interp)
{
    // fclaw_bool time_interp = alpha > 0; //
    int finer_level = level;
    int coarser_level = level - 1;

    exchange_info e_info;
    e_info.time_interp = time_interp;
    e_info.level = level;

    /* -----------------------------------------------------------
       Face average
       ----------------------------------------------------------- */

    // printf("Face average (mpirank = %d)\n",domain->mpirank);
    /* First pass : Iterate over coarse grids in space-filling curve */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 cb_face_average, (void *) &e_info);
    /* Second pass : Iterate over fine grids in space-filling curve,
       looking for parallel ghost patches or time interpolated data.
    */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain, finer_level,
                                 cb_face_average, (void *) &e_info);

    fclaw2d_domain_barrier(domain);

    /* -----------------------------------------------------------
       Corner average
       ----------------------------------------------------------- */

    /* First pass : Iterate over coarse grids in space filling curve */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;

    fclaw2d_domain_iterate_level(domain,coarser_level, cb_corner_average,
                                 (void *) &e_info);


    /* Second pass : Average over finer grids in sf curve */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain,finer_level, cb_corner_average,
                                 (void *) &e_info);

    /* -----------------------------------------------------------
       Physical BC : Set coarse grid physical boundary conditions -
       this will help with interpolation to finer grids. Time level
       't_level' is the time at the finer level, i.e.
       coarse_time + alpha*dt_coarse
       ----------------------------------------------------------- */
    set_phys_bc(domain,coarser_level,t_level,time_interp);

    /* -----------------------------------------------------------
       Face interpolate
       ----------------------------------------------------------- */
    /* First pass : Iterate over coarse grids in space filling curve */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain,coarser_level,cb_face_interpolate,
                                 (void *) &e_info);

    /* Second pass : Average over finer grids in sf curve */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;

    fclaw2d_domain_iterate_level(domain,finer_level,cb_face_interpolate,
                                 (void *) &e_info);
    /* -----------------------------------------------------------
       Corner interpolate
       ----------------------------------------------------------- */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain,coarser_level, cb_corner_interpolate,
                                 (void *) &e_info);

    /* Second pass : Iterate over fine grids on this processor that have
       off processor coarse grid neighbors.  The first pass doesn't catch
       these coarse grids, and the fine grids don't get ghost cell data set.
    */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain,finer_level, cb_corner_interpolate,
                                 (void *) &e_info);

}
