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

/* NOTE: Do we need the extern "C" here?  We're passing callbacks
   to C iterators but maybe C++ handles it just fine. */
#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

static
void cb_corner_average(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;
    bool intersects_bc[NumFaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    bool intersects_block[NumFaces];
    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    // Number of patch corners, not the number of corners in the domain!
    // const int numcorners = get_corners_per_patch(domain);

    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        // Get faces that intersect 'icorner'
        // There must be a clever way to do this...
        // p4est has tons of lookup table like this, can be exposed similarly
        int faces[SpaceDim];
        fclaw2d_domain_corner_faces(domain, icorner, faces);

        // Both faces are at a physical boundary
        bool is_phys_corner =
                intersects_bc[faces[0]] && intersects_bc[faces[1]];

        // Corner lies in interior of physical boundary edge.
        bool corner_on_phys_face = !is_phys_corner &&
                (intersects_bc[faces[0]] || intersects_bc[faces[1]]);

        bool interior_corner = !corner_on_phys_face && !is_phys_corner;

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
            bool is_block_corner =
                intersects_block[faces[0]] && intersects_block[faces[1]];

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
                // No corner neighbor
            }
            else if (ref_flag == 1)
            {
                fclaw2d_block_t *corner_block = &domain->blocks[corner_block_idx];
                fclaw2d_patch_t *corner_patch = &corner_block->patches[corner_patch_idx];
                ClawPatch *corner_cp = get_clawpatch(corner_patch);

                bool &time_interp = *((bool*) user);
                if (this_block_idx == corner_block_idx)
                {
                    this_cp->average_corner_ghost(icorner,refratio,corner_cp,time_interp);
                }
                else
                {
                    this_cp->mb_average_corner_ghost(icorner,refratio,corner_cp,time_interp,
                                                     is_block_corner,intersects_block);
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
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;
    bool intersects_bc[NumFaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    bool intersects_block[NumFaces];
    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    // Number of patch corners, not the number of corners in the domain!
    // const int numcorners = get_corners_per_patch(domain);

    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        // Get faces that intersect 'icorner'
        // There must be a clever way to do this...
        // p4est has tons of lookup table like this, can be exposed similarly
        int faces[SpaceDim];
        fclaw2d_domain_corner_faces(domain, icorner, faces);

        // Both faces are at a physical boundary
        bool is_phys_corner =
                intersects_bc[faces[0]] && intersects_bc[faces[1]];

        // Corner lies in interior of physical boundary edge.
        bool corner_on_phys_face = !is_phys_corner &&
                (intersects_bc[faces[0]] || intersects_bc[faces[1]]);

        bool interior_corner = !corner_on_phys_face && !is_phys_corner;

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
            bool is_block_corner =
                intersects_block[faces[0]] && intersects_block[faces[1]];

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
                // no corner neighbor
            }
            else if (ref_flag == 1)
            {
                // This can be cleaned up by just returning the neighbor patch directly, no?
                fclaw2d_block_t *corner_block = &domain->blocks[corner_block_idx];
                fclaw2d_patch_t *corner_patch = &corner_block->patches[corner_patch_idx];
                ClawPatch *corner_cp = get_clawpatch(corner_patch);

                bool &time_interp = *((bool*) user);
                if (this_block_idx == corner_block_idx)
                {
                    this_cp->interpolate_corner_ghost(icorner,refratio,corner_cp,time_interp);
                }
                else
                {
                    this_cp->mb_interpolate_corner_ghost(icorner,refratio,
                                                         corner_cp,time_interp,
                                                         is_block_corner, intersects_block);
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
    // const int p4est_refineFactor = get_p4est_refineFactor(domain);

    // We may need to average onto a time interpolated grid, not the actual solution.
    bool &time_interp = *((bool*) user);

    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    ClawPatch *this_cp = get_clawpatch(this_patch);

    for (int idir = 0; idir < 2; idir++)
    {
        // Loop over low side and high side
        for (int iface = 2*idir; iface <= 2*idir + 1; iface++)
        {
            int neighbor_block_idx;
            int neighbor_patch_idx[p4est_refineFactor];  // Be prepared to store 1 or more patch
                                                         // indices.
            int ref_flag;
            int *ref_flag_ptr = &ref_flag; // = -1, 0, 1

            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iface,
                               &neighbor_block_idx,
                               neighbor_patch_idx,
                               &ref_flag_ptr);

            if (ref_flag_ptr == NULL)
            {
                // no face neighor
            }
            else if (ref_flag == 1)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'this_patch' by averaging data from finer neighbors
                fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                ClawPatch *fine_neighbor_cp[p4est_refineFactor];
                for (int ir = 0; ir < p4est_refineFactor; ir++)
                {
                    fclaw2d_patch_t *neighbor_patch =
                        &neighbor_block->patches[neighbor_patch_idx[ir]];
                    fine_neighbor_cp[ir] = get_clawpatch(neighbor_patch);
                }
                bool block_boundary = this_block_idx != neighbor_block_idx;
                // Fill in ghost cells on 'this_cp' by averaging from 'fine_neighbor_cp'
                this_cp->average_face_ghost(idir,iface,p4est_refineFactor,refratio,
                                            fine_neighbor_cp,time_interp,block_boundary);
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}


// Iterator over patches looking for finer neighbors
static
void cb_face_interpolate(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    ClawPatch *this_cp = get_clawpatch(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iface = 2*idir; iface <= 2*idir + 1; iface++)
        {
            int neighbor_block_idx;
            int neighbor_patch_idx[p4est_refineFactor];
            int ref_flag;
            int *ref_flag_ptr = &ref_flag; // = -1, 0, 1

            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iface,
                               &neighbor_block_idx,
                               neighbor_patch_idx,
                               &ref_flag_ptr);

            if (ref_flag_ptr == NULL)
            {
                // no face neighbor - physical boundary?
            }
            else if (ref_flag == 1)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'neighbor_patch' by interpolation
                fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                ClawPatch *fine_neighbor_cp[p4est_refineFactor];
                for (int ir = 0; ir < p4est_refineFactor; ir++)
                {
                    fclaw2d_patch_t *neighbor_patch =
                        &neighbor_block->patches[neighbor_patch_idx[ir]];
                    fine_neighbor_cp[ir] = get_clawpatch(neighbor_patch);
                }

                // Fill in fine grid ghost on 'fine_neighbor_cp' by interpolating from 'this_cp',
                // doing time interpolation if necessary
                bool &time_interp = *((bool*) user);
                bool block_boundary = this_block_idx != neighbor_block_idx;
                this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                refratio,fine_neighbor_cp,time_interp,
                                                block_boundary);
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}


static
void cb_setup_time_interp(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    // This is called for all patches at a level coarser than the one we have.
    ClawPatch *cp = get_clawpatch(this_patch);

    double &alpha = *((double*) user);
    cp->time_interpolate(alpha);
}

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

// ----------------------------------------------------------------------
// Main routine in this file
// ----------------------------------------------------------------------

void exchange_with_coarse(fclaw2d_domain_t *domain,
                          int level, double t_level,
                          double alpha)
{
    // Simple exchange - no time interpolation needed
    bool time_interp = alpha > 0; //
    int coarser_level = level - 1;

    if (time_interp)
    {
        fclaw2d_domain_iterate_level(domain, coarser_level,
                                     cb_setup_time_interp,
                                     (void *) &alpha);
    }


/* -------------------------------------------------------------------
   Fill coarse grid ghost cells at edges and corners that are shared with
   the finer grid
   -----------------------------------------------------------------------*/

    // Iterate over coarser level and average from finer neighbors to coarse.
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 cb_face_average,
                                 (void *) &time_interp);

    // Average fine grid corners to the coarse grid ghost cells
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 cb_corner_average,
                                 (void *) &time_interp);

    // Set coarse grid physical boundary conditions - this will help with
    // interpolation to finer grids.
    set_phys_bc(domain,coarser_level,t_level);

    // Interpolate coarse grid to fine.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 cb_face_interpolate,
                                 (void *) &time_interp);

    // Interpolate coarse grid to fine grid ghost cells.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 cb_corner_interpolate,
                                 (void *) &time_interp);
}
