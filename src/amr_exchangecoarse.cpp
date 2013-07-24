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

            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                &corner_block_idx,
                                &neighbor_patch,
                                &ref_flag_ptr,
                                is_block_corner);

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
                    int icorner_coarse = 3-icorner;

                    if (this_block_idx == corner_block_idx)
                    {
                        coarse_cp->average_corner_ghost(icorner_coarse,
                                                      refratio,corner_cp,time_interp);
                    }
                    else
                    {
                        coarse_cp->mb_average_corner_ghost(icorner_coarse,
                                                           refratio,corner_cp,time_interp,
                                                           is_block_corner,intersects_block);
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

            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                &corner_block_idx,
                                &neighbor_patch,
                                &ref_flag_ptr,
                                is_block_corner);

            if (ref_flag_ptr == NULL)
            {
                // no corner neighbor
            }
            else if (ref_flag == 1)
            {
                ClawPatch *corner_cp = get_clawpatch(neighbor_patch);

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
                    int icorner_coarse = 3-icorner;

                    if (this_block_idx == corner_block_idx)
                    {
                        coarse_cp->interpolate_corner_ghost(icorner_coarse,
                                                            refratio,corner_cp,time_interp);
                    }
                    else
                    {
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

    ClawPatch *this_cp = get_clawpatch(this_patch);

    for (int idir = 0; idir < 2; idir++)
    {
        // Loop over low side and high side
        for (int iface = 2*idir; iface <= 2*idir + 1; iface++)
        {
            int neighbor_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag; // = -1, 0, 1
            fclaw2d_patch_t *neighbor_patches[p4est_refineFactor];

            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iface,
                               &neighbor_block_idx,
                               neighbor_patches,
                               &ref_flag_ptr);

            if (ref_flag_ptr == NULL)
            {
                // no face neighor
            }
            else if (ref_flag == 1 && is_coarse)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'this_patch' by averaging data from finer neighbors
                fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
                for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                {
                    ClawPatch* fine_neighbor_cp = get_clawpatch(neighbor_patches[igrid]);
                    this_cp->average_face_ghost(idir,iface,p4est_refineFactor,refratio,
                                                fine_neighbor_cp,time_interp,block_boundary,
                                                igrid);
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
                    int igrid;
                    if (idir == 0)
                    {
                        /* this is a awkward;  is there a better way to do this? */
                        double tol = 1e-12;
                        double ylow_fine = fine_cp->ylower();
                        double yhi_fine = fine_cp->yupper();
                        double ylow_coarse = coarse_cp->ylower();
                        double yhi_coarse = coarse_cp->yupper();
                        if (fabs(ylow_fine - ylow_coarse) < tol)
                            igrid = 0;
                        else if (fabs(yhi_fine - yhi_coarse) < tol)
                            igrid = 1;
                        else
                        {
                            printf("amr_exchange_coarse; cannot find igrid\n");
                            exit(1);
                        }
                    }
                    else
                    {
                        /* this is a awkward;  is there a better way to do this? */
                        double tol = 1e-12;
                        double xlow_fine = fine_cp->xlower();
                        double xhi_fine = fine_cp->xupper();
                        double xlow_coarse = coarse_cp->xlower();
                        double xhi_coarse = coarse_cp->xupper();
                        if (fabs(xlow_fine - xlow_coarse) < tol)
                            igrid = 0;
                        else if (fabs(xhi_fine - xhi_coarse) < tol)
                            igrid = 1;
                        else
                        {
                            printf("amr_exchange_coarse; cannot find igrid\n");
                            exit(1);
                        }
                    }
                    /* Swap out coarse and fine */
                    int iface_coarse;
                    if (iface == 0)
                        iface_coarse = 1;
                    else if (iface == 1)
                        iface_coarse = 0;
                    else if (iface == 2)
                        iface_coarse = 3;
                    else if (iface == 3)
                        iface_coarse = 2;

                    fclaw_bool block_boundary = this_block_idx != neighbor_block_idx;
                    coarse_cp->average_face_ghost(idir,iface_coarse,p4est_refineFactor,refratio,
                                                fine_cp,time_interp,block_boundary,
                                                igrid);

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

    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    ClawPatch *this_cp = get_clawpatch(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iface = 2*idir; iface <= 2*idir + 1; iface++)
        {
            int neighbor_block_idx;
            int ref_flag;
            int *ref_flag_ptr = &ref_flag; // = -1, 0, 1
            fclaw2d_patch_t *neighbor_patches[p4est_refineFactor];

            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iface,
                               &neighbor_block_idx,
                               neighbor_patches,
                               &ref_flag_ptr);

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
                    this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                    refratio,fine_neighbor_cp,time_interp,
                                                    block_boundary,igrid);
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
                    int igrid;
                    if (idir == 0)
                    {
                        /* this is a awkward;  is there a better way to do this? */
                        double tol = 1e-12;
                        double ylow_fine = fine_cp->ylower();
                        double yhi_fine = fine_cp->yupper();
                        double ylow_coarse = coarse_cp->ylower();
                        double yhi_coarse = coarse_cp->yupper();
                        if (fabs(ylow_fine - ylow_coarse) < tol)
                            igrid = 0;
                        else if (fabs(yhi_fine - yhi_coarse) < tol)
                            igrid = 1;
                        else
                        {
                            printf("amr_exchange_coarse; cannot find igrid\n");
                            exit(1);
                        }
                    }
                    else
                    {
                        /* this is a awkward;  is there a better way to do this? */
                        double tol = 1e-12;
                        double xlow_fine = this_cp->xlower();
                        double xhi_fine = this_cp->xupper();
                        double xlow_coarse = coarse_cp->xlower();
                        double xhi_coarse = coarse_cp->xupper();
                        if (fabs(xlow_fine - xlow_coarse) < tol)
                            igrid = 0;
                        else if (fabs(xhi_fine - xhi_coarse) < tol)
                            igrid = 1;
                        else
                        {
                            printf("amr_exchange_coarse; cannot find igrid\n");
                            exit(1);
                        }
                    }
                    /* Swap out coarse and fine */
                    int iface_coarse;
                    if (iface == 0)
                        iface_coarse = 1;
                    else if (iface == 1)
                        iface_coarse = 0;
                    else if (iface == 2)
                        iface_coarse = 3;
                    else if (iface == 3)
                        iface_coarse = 2;

                    coarse_cp->interpolate_face_ghost(idir,iface_coarse,p4est_refineFactor,refratio,
                                                      fine_cp,time_interp,block_boundary,
                                                      igrid);

                }
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
    /* This is called for all patches on the coarse level */
    ClawPatch *cp = get_clawpatch(this_patch);
    double &alpha = *((double*) user);

    /* This constructs a time interpolated version of the data on
       the coarser grid */
    cp->setup_for_time_interpolation(alpha);
}



/* ----------------------------------------------------------------------
   Main routine in this file.  This file assumes that both coarse and
   fine grids have valid interior data;  they only need to exchange (
   via interpolating and averaging) ghost cell values.
   -------------------------------------------------------------------- */

void exchange_with_coarse(fclaw2d_domain_t *domain,
                          int level, double t_level,
                          double alpha)
{
    exchange_info e_info;
    // Simple exchange - no time interpolation needed
    fclaw_bool time_interp = alpha > 0; //
    int finer_level = level;
    int coarser_level = level - 1;

    e_info.time_interp = time_interp;
    e_info.level = level;

    /* copy griddata to griddata_time_sync so that parallel exchange will work */
    if (time_interp)
    {
        fclaw2d_domain_iterate_level(domain, coarser_level,cb_setup_time_interp,
                                     (void *) &alpha);
    }

    /* Do parallel ghost exchange; Pointers to data get reassigned here, in case
       we are in the time interpolated case */
    exchange_ghost_patch_data(domain,time_interp);


    /* Iterate over coarser level and average from finer neighbors to coarse. */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 cb_face_average, (void *) &e_info);

    /* Iterate over coarser level and average from finer neighbors to coarse. */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain, finer_level,
                                 cb_face_average, (void *) &e_info);

    /* Average fine grid corners to the coarse grid ghost cells.  This will pick
     any finer grid ghost patches and average them onto on-proc coarse grids*/
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain,coarser_level, cb_corner_average,
                                 (void *) &e_info);

    /* Average fine grid corners to the coarse grid ghost cells.  This will pick
     any finer grid ghost patches and average them onto on-proc coarse grids*/
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain,finer_level, cb_corner_average,
                                 (void *) &e_info);


    /* Set coarse grid physical boundary conditions - this will help with
       interpolation to finer grids. Time level 't_level' is the time
       at the finer level, i.e. coarse_time + alpha*dt_coarse
    */
    set_phys_bc(domain,coarser_level,t_level,time_interp);

    /* Interpolate from coarse grids to finer grids. */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain,coarser_level,cb_face_interpolate,
                                 (void *) &e_info);

    /* Interpolate from coarse grids to finer grids. */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain,finer_level,cb_face_interpolate,
                                 (void *) &e_info);

    /* Interpolate coarse grid to fine grid ghost cells. */
    e_info.is_coarse = fclaw_true;
    e_info.is_fine = fclaw_false;
    fclaw2d_domain_iterate_level(domain,coarser_level, cb_corner_interpolate,
                                 (void *) &e_info);

    /* Interpolate coarse grid to fine grid ghost cells. */
    e_info.is_coarse = fclaw_false;
    e_info.is_fine = fclaw_true;
    fclaw2d_domain_iterate_level(domain,finer_level, cb_corner_interpolate,
                                 (void *) &e_info);

}
