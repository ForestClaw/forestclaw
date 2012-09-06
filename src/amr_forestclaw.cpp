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
#include "fclaw2d_convenience.h"

#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>

class ClawPatch;
class global_parms;

// -----------------------------------------------------------------
// Setting data in domain and patches
// -----------------------------------------------------------------

/* Proposed naming convention:
 * Functions called *init* allocate memory inside an object.
 * Functions called *new* create a new object and then does an init.
 * Functions called *set* update a memory location without allocating.
 * Functions called *get* return a memory location without (de)allocating.
 * Functions called *reset* deallocate memory in a given object.
 * Functions called *destroy* do a reset and then destroy the object itself.
 */

/* domain data */

static void init_domain_data(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = FCLAW2D_ALLOC_ZERO (fclaw2d_domain_data_t, 1);
    domain->user = (void *) ddata;
}

void set_domain_data(fclaw2d_domain_t *domain, global_parms *gparms, const amr_options_t *amropts)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    ddata->parms = gparms;
    ddata->amropts = amropts;
}


fclaw2d_domain_data_t *get_domain_data(fclaw2d_domain_t *domain)
{
    return (fclaw2d_domain_data_t *) domain->user;
}

/* block data */

static void init_block_data(fclaw2d_block_t *block)
{
    fclaw2d_block_data_t *bdata = FCLAW2D_ALLOC_ZERO (fclaw2d_block_data_t, 1);
    block->user = (void *) bdata;
}

void set_block_data(fclaw2d_block_t *block, int *mthbc)
{
    fclaw2d_block_data_t *bdata = get_block_data(block);
    for(int i = 0; i < NumFaces; i++)
    {
        bdata->mthbc[i] = mthbc[i];
    }
}


fclaw2d_block_data_t *get_block_data(fclaw2d_block_t *block)
{
    return (fclaw2d_block_data_t *) block->user;
}

/* patch data */

void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp)
{
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;

    pdata->cp = cp;
}


/* These are lazy helper functions that really do make things easier.. */
global_parms* get_domain_parms(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    return ddata->parms;
}

void set_domain_time(fclaw2d_domain_t *domain, Real time)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    ddata->curr_time = time;
}

Real get_domain_time(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    return ddata->curr_time;
}


// Will change the name of this to 'get_clawpatch' eventually
ClawPatch* get_patch_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata;
    pdata = (fclaw2d_patch_data_t *) patch->user;

    return pdata->cp;
}
/* end of helper functions */



void set_problem_parameters()
{
    setprob_();
}



// -----------------------------------------------------------------
// Diagnostics
// -----------------------------------------------------------------
void cb_check_conservation(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           void *user)
{
    ClawPatch *this_cp = get_patch_data(this_patch);
    Real *sum = (Real*) user;

    *sum += this_cp->compute_sum();
}


void check_conservation(fclaw2d_domain_t *domain)
{
    Real sum = 0;
    fclaw2d_domain_iterate_patches(domain,cb_check_conservation,(void *) &sum);

    printf("Total sum = %24.16f\n",sum);
}




// -----------------------------------------------------------------
// Get face and corner neighbors
// -----------------------------------------------------------------
void get_face_neighbors(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        int iside,
                        int *neighbor_block_idx,
                        int neighbor_patch_idx[],
                        int *ref_flag,
                        bool *is_phys_bc)
{
    const int num_face_corners = fclaw2d_domain_num_face_corners (domain);
    int rproc[num_face_corners];
    int rblockno;
    int rpatchno[num_face_corners];
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
        *ref_flag = 0;  // Want to have a valid value for 'ref_flag', so that it can be checked
                        // outside of this routine. Physical BCs ghost cells are in fact at the
                        // same  level as 'this_patch', so seems reasonable to set it to 0.
        *is_phys_bc = true;
     }
    else
    {
        *is_phys_bc = false;
        if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            // Neighbors are finer grids
            *ref_flag = 1;
            for(int ir = 0; ir < num_face_corners; ir++)
            {
                neighbor_patch_idx[ir] = rpatchno[ir];
            }
        }
        else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            // Neighbor is at the same level
            *ref_flag = 0;
            neighbor_patch_idx[0] = rpatchno[0];
        }
        else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
        {
            // Neighbor is a coarser grid
            *ref_flag = -1;
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
                         int *ref_flag)
{
    int neighbor_corner = fclaw2d_patch_corner_neighbors(domain, this_block_idx, this_patch_idx, icorner);

    if (neighbor_corner == -1)
    {
        printf("Patch %d at corner %d does not have any corner neighbors\n",
               this_patch_idx,icorner);
        exit(1);
    }

    *corner_block_idx = this_block_idx;
    *corner_patch_idx = neighbor_corner;
    *ref_flag = 0; // only return patches at the same level for now
}

void get_phys_boundary(fclaw2d_domain_t *domain,
                       int this_block_idx,
                       int this_patch_idx,
                       bool *intersects_bc)
{
    int bdry[NumFaces];
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);

    for(int i = 0; i < NumFaces; i++)
    {
        intersects_bc[i] = bdry[i] > 0;
    }
}


// -----------------------------------------------------------------
// Exchange corner and face information
// -----------------------------------------------------------------
void cb_bc_level_face_exchange(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               void *user)
{
    const int num_face_corners = fclaw2d_domain_num_face_corners (domain);

    ClawPatch *this_cp = get_patch_data(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        int iface = 2*idir + 1;  // Look only at high side for level copy.

        // Output arguments
        int neighbor_block_idx;
        int neighbor_patch_idx[num_face_corners];  // Be prepared to store 1 or more patch indices.
        int ref_flag; // = -1, 0, 1
        bool is_phys_bc;
        get_face_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iface,
                           &neighbor_block_idx,
                           neighbor_patch_idx,
                           &ref_flag,
                           &is_phys_bc);


        if (ref_flag == 0 && !is_phys_bc)
        {
            // Ghost cells overlap patch at the same level, and so we can do a straight copy.
            fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
            fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[0]];
            ClawPatch *neighbor_cp = get_patch_data(neighbor_patch);

            // Exchange between 'this_patch' and 'neighbor patch(es)' in direction 'idir'
            this_cp->exchange_face_ghost(idir,neighbor_cp);
        }
    } // loop over directions (idir = 0,1,2)


    // Now check for any physical boundary conditions on this patch
    bool intersects_bc[NumFaces];
    Real curr_time = get_domain_time(domain);
    Real dt = 1e20;   // When do we need dt in setting a boundary condition?
    get_phys_boundary(domain,this_block_idx,this_patch_idx,intersects_bc);


    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    fclaw2d_block_data_t *bdata = get_block_data (this_block);
    this_cp->set_phys_face_ghost(intersects_bc,bdata->mthbc,curr_time,dt);

    // Do this twice to get corners correctly.
    this_cp->set_phys_face_ghost(intersects_bc,bdata->mthbc,curr_time,dt);
}

void cb_level_corner_exchange(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              void *user)
{
    const int num_face_corners = fclaw2d_domain_num_face_corners (domain);
    bool intersects_bc[NumFaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

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

        ClawPatch *this_cp = get_patch_data(this_patch);
        if (is_phys_corner)
        {
            // This corner has no neighbor patches;
            // try to do something sensible using only what
            // we know about this physical boundary conditions at this corner
            fclaw2d_block_t *block = &domain->blocks[this_block_idx];
            fclaw2d_block_data_t *bdata = get_block_data (block);
            Real curr_time = get_domain_time(domain);
            Real dt = 1e20;   // When do we need dt in setting a boundary condition?
            this_cp->set_phys_corner_ghost(icorner,bdata->mthbc,curr_time,dt);
        }
        else if (corner_on_phys_face)
        {
            // Corners lies some where along a physical edge,
            // but not at a physical corner.
            int iside;
            if (intersects_bc[faces[0]])
            {
                iside = faces[1];
            }
            else
            {
                iside = faces[0];
            }
            if (iside % 2 == 1) // hi-side face
            {
                // only initiate exchange with a hi-side neighbor
                // int refratio = gparms->m_refratio;
                int neighbor_block_idx;
                int neighbor_patch_idx[num_face_corners];  // Be prepared to store 1 or more patch indices.
                int ref_flag; // = -1, 0, 1
                bool is_phys_bc;
                get_face_neighbors(domain,
                                   this_block_idx,
                                   this_patch_idx,
                                   iside,
                                   &neighbor_block_idx,
                                   neighbor_patch_idx,
                                   &ref_flag,
                                   &is_phys_bc);
                if (ref_flag == 0)
                {
                    if (is_phys_bc)
                    {
                        fprintf (stderr, "ERROR: For phys bc faces patch_idx is undefined\n");
                        exit (1);
                    }

                    // Only doing a level exchange now
                    fclaw2d_block_t *neighbor_block =
                        &domain->blocks[neighbor_block_idx];
                    fclaw2d_patch_t *neighbor_patch =
                        &neighbor_block->patches[neighbor_patch_idx[0]];
                    ClawPatch *neighbor_cp = get_patch_data(neighbor_patch);
                    this_cp->exchange_phys_face_corner_ghost(icorner,iside,
                        neighbor_cp);
                }
            }
        }
        else if (interior_corner)
        {
            // The hi-side faces above are only those faces with an endpoint on the
            // physical boundary, i.e. which intersect a corner that lies on an physical
            // boundary.  Each one of these corners is on the 'hi-side' of another patch
            // face.
            //
            // Only initiate exchange if we are at a high side corner
            if (icorner % 2 == 1) // only exchange with high side corners (corners 1, 3, 5 or 7)
            {
                int corner_block_idx;
                int corner_patch_idx;
                int ref_flag;


                get_corner_neighbor(domain,
                                    this_block_idx,
                                    this_patch_idx,
                                    icorner,
                                    &corner_block_idx,
                                    &corner_patch_idx,
                                    &ref_flag);

                if (ref_flag == 0)
                {
                    // Upper right of 'this_cp' exchanges with lower left of 'corner_cp' or
                    // lower right of 'this_cp' exchanges with upper left of 'corner_cp'

                    // Only consider case in which neighbor is at same level of refinement
                    fclaw2d_block_t *corner_block = &domain->blocks[corner_block_idx];
                    fclaw2d_patch_t *corner_patch = &corner_block->patches[corner_patch_idx];
                    ClawPatch *corner_cp = get_patch_data(corner_patch);

                    this_cp->exchange_corner_ghost(icorner,corner_cp);
                }
            }
        }
    }
}


void bc_level_exchange(fclaw2d_domain_t *domain, int a_level)
{

    fclaw2d_domain_iterate_level(domain, a_level,
                                  cb_bc_level_face_exchange, (void *) NULL);

    // Do corner exchange only after physical boundary conditions have been set on all patches,
    // since corners may overlap phyical ghost cell region of neighboring patch.
    fclaw2d_domain_iterate_level(domain, a_level, cb_level_corner_exchange, (void *) NULL);
}


// -----------------------------------------------------------------
// Multi-level ghost cell operations
//   -- average fine grid ghost cells to coarse grid
//   -- interpolate coarse grid to fine grid ghost cells
// -----------------------------------------------------------------
void cb_bc_average(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    const int num_face_corners = fclaw2d_domain_num_face_corners (domain);

    // We may need to average onto a time interpolated grid, not the actual solution.
    fclaw2d_subcycle_info *step_info = (fclaw2d_subcycle_info*) user;
    bool do_time_interp = step_info->do_time_interp;

    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    global_parms *gparms = ddata->parms;
    int refratio = gparms->m_refratio;

    ClawPatch *this_cp = get_patch_data(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iside = 2*idir; iside < 2*idir + 1; iside++)
        {
            int neighbor_block_idx;
            int neighbor_patch_idx[num_face_corners];  // Be prepared to store 1 or more patch indices.
            int ref_flag; // = -1, 0, 1
            bool is_phys_bc;

            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iside,
                               &neighbor_block_idx,
                               neighbor_patch_idx,
                               &ref_flag,
                               &is_phys_bc);

            if (ref_flag == 1)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'this_patch' by averaging data from finer neighbors
                fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                ClawPatch *fine_neighbor_cp[num_face_corners];
                for (int ir = 0; ir < num_face_corners; ir++)
                {
                    fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[ir]];
                    fine_neighbor_cp[ir] = get_patch_data(neighbor_patch);
                }

                // Fill in ghost cells on 'this_cp' by averaging from 'fine_neighbor_cp'
                this_cp->average_face_ghost(idir,iside,num_face_corners,
                                            refratio,fine_neighbor_cp,do_time_interp);
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}


void cb_bc_interpolate(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    const int num_face_corners = fclaw2d_domain_num_face_corners (domain);
    fclaw2d_subcycle_info_t *step_info = (fclaw2d_subcycle_info_t*) user;

    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    global_parms *gparms = ddata->parms;
    int refratio = gparms->m_refratio;

    ClawPatch *this_cp = get_patch_data(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iside = 2*idir; iside < 2*idir + 1; iside++)
        {
            int neighbor_block_idx;
            int neighbor_patch_idx[num_face_corners];
            int ref_flag; // = -1, 0, 1
            bool is_phys_bc;

            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iside,
                               &neighbor_block_idx,
                               neighbor_patch_idx,
                               &ref_flag,
                               &is_phys_bc);

            if (ref_flag == 1)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'neighbor_patch' by interpolation
                fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                ClawPatch *fine_neighbor_cp[num_face_corners];
                for (int ir = 0; ir < num_face_corners; ir++)
                {
                    fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[ir]];
                    fine_neighbor_cp[ir] = get_patch_data(neighbor_patch);
                }

                // Fill in fine grid ghost on 'fine_neighbor_cp' by interpolating from 'this_cp',
                // doing time interpolation if necessary
                this_cp->interpolate_face_ghost(idir,iside,num_face_corners,
                                                refratio,fine_neighbor_cp,step_info->do_time_interp);
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}

void cb_setup_time_interp(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    // construct all coarse level time interpolated intermediate grids.  Interpolate ghost
    // values as well, even though neighboring fine grids may overwrite some ghost values
    fclaw2d_subcycle_info_t *step_info = (fclaw2d_subcycle_info_t*) user;

    ClawPatch *cp = get_patch_data(this_patch);
    cp->time_interpolate(step_info->fine_step, step_info->coarse_step, step_info->refratio);
}

void bc_exchange_with_coarse_time_interp(fclaw2d_domain_t *domain, const int& a_level, const int& a_coarse_step,
                                         const int& a_fine_step,const int& a_refratio)
{
    // First, average fine grid to coarse grid cells
    fclaw2d_subcycle_info_t step_info;
    step_info.coarse_step = a_coarse_step;
    step_info.fine_step = a_fine_step;
    step_info.refratio = a_refratio;
    step_info.do_time_interp = true;

    // Set up patch data for time interpolation.
    int coarser_level = a_level - 1;
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_setup_time_interp,
                                 (void *) &step_info);


    // Average onto time interpolated 'virtual' data so that we can use these averaged ghost
    // values in the interpolation step below.
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_average,
                                 (void *) &step_info);

    // Interpolate coarse grid to fine.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_interpolate,
                                 (void *) &step_info);
}

void bc_exchange_with_coarse(fclaw2d_domain_t *domain, const int& a_level)
{
    // Simple exchange - no time interpolation needed
    fclaw2d_subcycle_info_t step_info;
    step_info.do_time_interp = false;

    int coarser_level = a_level - 1;
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_average,
                                 (void *) &step_info);

    // Interpolate coarse grid to fine.  This may require time interpolation
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_interpolate,
                                 (void *) &step_info);
}


// -----------------------------------------------------------------
// Time stepping
//   -- saving time steps
//   -- restoring time steps
//   -- advancing levels
// -----------------------------------------------------------------
void cb_restore_time_step(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    ClawPatch *this_cp = get_patch_data(this_patch);

    // Copy most current time step data to grid data. (m_griddata <== m_griddata_last)
    this_cp->restore_step();
}


void restore_time_step(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_iterate_patches(domain,cb_restore_time_step,(void *) NULL);
}

void cb_save_time_step(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    ClawPatch *this_cp = get_patch_data(this_patch);

    // Copy grid data (m_griddata) on each patch to temporary storage (m_griddata_tmp <== m_griddata);
    this_cp->save_step();
}


void save_time_step(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_iterate_patches(domain,cb_save_time_step,(void *) NULL);
}


void cb_advance_patch(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    global_parms *gparms = ddata->parms;

    ClawPatch *cp = get_patch_data(this_patch);
    fclaw2d_level_time_data_t *time_data = (fclaw2d_level_time_data_t *) user;

    Real dt = time_data->dt;
    Real t = time_data->t;

    int level = this_patch->level;

    Real maxcfl_grid = cp->step_noqad(t,dt,level,*gparms);

    time_data->maxcfl = max(maxcfl_grid,time_data->maxcfl);
}


Real advance_level(fclaw2d_domain_t *domain,
                   const int& a_level,
                   const int& a_curr_fine_step,
                   subcycle_manager* a_time_stepper)
{

    // printf("Advancing level %d by dt = %16.4e\n",a_level,a_time_stepper.dt(a_level));
    // Check BCs
    Real maxcfl = 0;
    if (!a_time_stepper->can_advance(a_level,a_curr_fine_step))
    {
        if (!a_time_stepper->level_exchange_done(a_level))
        {
            // Level exchange should have been done right after solution update.
            printf("Error (advance_level) : Level exchange at level %d not done at time step %d\n",a_level,a_curr_fine_step);
            exit(1);
        }
        if (!a_time_stepper->exchanged_with_coarser(a_level))
        {
            printf("For uniform refinement, I should not be here!\n");
            exit(1);
            int last_coarse_step = a_time_stepper->last_step(a_level-1);
            if (a_curr_fine_step == last_coarse_step)
            {
                // Levels are time synchronized and we can do a simple coarse/fine
                // exchange without time interpolation or advancing the coarser level
                bc_exchange_with_coarse(domain,a_level);
                a_time_stepper->increment_coarse_exchange_counter(a_level);
                a_time_stepper->increment_fine_exchange_counter(a_level-1);
            }
            else
            {
                if ((a_curr_fine_step > last_coarse_step) ||
                    (a_time_stepper->nosubcycle() && !a_time_stepper->is_coarsest(a_level)))
                {
                    // (1) subcycled case : a_curr_fine_step will only be greater than
                    // last_coarse_step if we haven't yet advanced the coarse grid to a time level
                    // beyond the current fine grid level.
                    //
                    // (2) non-subcycled case : this advance is a bit gratuitous, because we
                    // don't need it to advance the fine grid;  rather, we put this advance
                    // here as a way to get advances on the coarser levels.
                    maxcfl = advance_level(domain,a_level-1,last_coarse_step,a_time_stepper);
                }
                if (!a_time_stepper->nosubcycle())
                {
                    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
                    global_parms *gparms = ddata->parms;
                    int refratio = gparms->m_refratio;

                    // (1) a_curr_fine_step > last_coarse_step : we just advanced the coarse grid
                    // and can now apply time interpolated boundary conditions, or
                    //
                    // (2) a_curr_fine_step < last_coarse_step : we advanced the coarse
                    // grid in a previous step but we still have to do time interpolation (this
                    // can only happen if refratio > 2)

                    bc_exchange_with_coarse_time_interp(domain,a_level,last_coarse_step,a_curr_fine_step,refratio);
                    a_time_stepper->increment_coarse_exchange_counter(a_level);

                    // Don't increment the fine_exchange_counter, since in the time interpolated case,
                    // there is no coarse time data at time step a_curr_fine_step.
                    // a_time_stepper->increment_fine_exchange_counter(a_level-1);
                }
            }
        }
    }

    fclaw2d_level_time_data_t time_data;

    time_data.maxcfl = maxcfl;
    time_data.dt = a_time_stepper->dt(a_level);
    time_data.t = a_time_stepper->current_time(a_level);

    // Advance this level from 'a_curr_fine_step' to 'a_curr_fine_step + a_time_stepper.step_inc(a_level)'
    fclaw2d_domain_iterate_level(domain, a_level,
                                 (fclaw2d_patch_callback_t) cb_advance_patch,
                                 (void *) &time_data);

    a_time_stepper->increment_step_counter(a_level);
    a_time_stepper->increment_time(a_level);

    bc_level_exchange(domain,a_level);
    a_time_stepper->increment_level_exchange_counter(a_level);

    return time_data.maxcfl;  // Maximum from level iteration

}


Real advance_all_levels(fclaw2d_domain_t *domain,
                        subcycle_manager *a_time_stepper)
{
    // 'n_fine_steps' is the number of steps we must take on the finest level to equal one
    // step on the coarsest non-empty level, i.e. minlevel.
    int minlevel = a_time_stepper->minlevel();
    int n_fine_steps = a_time_stepper->step_inc(minlevel); // equal 1 in the non-subcycled case.
    int maxlevel = a_time_stepper->maxlevel();
    Real maxcfl = 0;
    for(int nf = 0; nf < n_fine_steps; nf++)
    {
        Real cfl_step = advance_level(domain,maxlevel,nf,a_time_stepper);
        maxcfl = max(cfl_step,maxcfl);
    }
    return maxcfl;
}


// -----------------------------------------------------------------
// Regridding
//   -- Initialization routines
//   -- cell tagging
//   -- interpolating/coarsening as needed
// -----------------------------------------------------------------

void cb_tag_patch(fclaw2d_domain_t *domain,
                  fclaw2d_patch_t *this_patch,
                  int this_block_idx,
                  int this_patch_idx,
                  void *user)
{
    ClawPatch *cp = get_patch_data(this_patch);

    bool level_refined = *((bool *) user);
    bool patch_refined = cp->tag_for_refinement();

    if (patch_refined)
    {
        fclaw2d_patch_mark_refine(domain, this_block_idx, this_patch_idx);
    }

    (*(bool *) user) =  level_refined || patch_refined;
}


void amr_set_base_level(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (domain);
    global_parms *gparms = ddata->parms;
    gparms->print_inputParams();

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        init_block_data(block);
        set_block_data(block,gparms->m_mthbc);

        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = new ClawPatch();

            cp->define(patch->xlower,
                       patch->ylower,
                       patch->xupper,
                       patch->yupper,
                       gparms);
            set_patch_data(patch,cp);
        }
    }
}


void cb_domain_adapt(fclaw2d_domain_t * old_domain,
                     fclaw2d_patch_t * old_patch,
                     fclaw2d_domain_t * new_domain,
                     fclaw2d_patch_t * new_patch,
                     fclaw2d_patch_relation_t newsize,
                     int blockno,
                     int old_patchno, int new_patchno,
                     void *user)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(old_domain);
    global_parms *gparms = ddata->parms;

    const int num_siblings = fclaw2d_domain_num_corners (old_domain);
    bool init_grid = *(bool *) user;

    int refratio = gparms->m_refratio;
    int maxlevel = gparms->m_maxlevel;

    if (newsize == FCLAW2D_PATCH_SAMESIZE)
    {
        // Grid was not coarsened or refined.
        ClawPatch *cp_old = get_patch_data(&old_patch[0]);
        ClawPatch *cp_new = new ClawPatch();
        cp_new->copy(cp_old);  // Copy grid data and aux data
        set_patch_data(&new_patch[0],cp_new);
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        // New grids are FINER grids
        ClawPatch *cp_old = get_patch_data(&old_patch[0]);

        for (int igrid = 0; igrid < num_siblings; igrid++)
        {
            ClawPatch *cp_new = new ClawPatch();

            cp_new->define(new_patch[igrid].xlower,
                           new_patch[igrid].ylower,
                           new_patch[igrid].xupper,
                           new_patch[igrid].yupper,
                           gparms);

            if (init_grid)
            {
                if (gparms->m_maux > 0)
                {
                    cp_new->setAuxArray(maxlevel,refratio,new_patch[igrid].level);
                }
                cp_new->initialize();
            }
            else
            {
                cp_old->interpolate_to_fine_patch(cp_new,igrid,p4est_refineFactor,refratio);
            }
            set_patch_data(&new_patch[igrid],cp_new);
        }
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        // newsize == DOUBLESIZE (must remember  : DOUBLESIZE means a coarser grid!)
        // new grid is a COARSE grid
        ClawPatch *cp_new = new ClawPatch();
        cp_new->define(new_patch[0].xlower,
                       new_patch[0].ylower,
                       new_patch[0].xupper,
                       new_patch[0].yupper,
                       gparms);
        set_patch_data(&new_patch[0],cp_new);

        for(int igrid = 0; igrid < num_siblings; igrid++)
        {
            // Get one of the older finer grids and average to new coarse grid.
            // Assume that we will never coarsen when we are initializing the grids.
            ClawPatch *cp_old = get_patch_data(&old_patch[igrid]);
            cp_new->coarsen_from_fine_patch(cp_old, igrid, p4est_refineFactor,refratio);
        }
    }
    else
    {
        printf("cb_adapt_domain : newsize not recognized\n");
        exit(1);
    }
}



// -----------------------------------------------------------------
// Initial grid
// -----------------------------------------------------------------
void cb_amrinit(fclaw2d_domain_t *domain,fclaw2d_patch_t *this_patch,
                int this_block_idx, int this_patch_idx, void *user)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    global_parms *gparms = ddata->parms;
    int refratio = gparms->m_refratio;
    int maxlevel = gparms->m_maxlevel;
    int maux = gparms->m_maux;

    ClawPatch *cp = get_patch_data(this_patch);

    if (maux > 0)
    {
        cp->setAuxArray(maxlevel,refratio,this_patch->level);
    }
    cp->initialize();
}

// Initialize a base level of grids
void amrinit(fclaw2d_domain_t **domain,
             global_parms *gparms, const amr_options_t * amropts)
{
    Real t = 0;
    cout << "domain (amrinit,beginning) " << *domain << endl;

    init_domain_data(*domain);
    set_domain_data(*domain, gparms, amropts);
    set_domain_time(*domain,t);

    // fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    // global_parms *gparms = ddata->parms;
    int minlevel = gparms->m_minlevel;
    int maxlevel = gparms->m_maxlevel;

    // Set problem dependent parameters for Riemann solvers, etc.
    // Values are typically stored in Fortran common blocks, and are not
    // available outside of Fortran.
    set_problem_parameters();

    // Set up storage for base level grids so we can initialize them
    // Allocates per-block and per-patch user data
    amr_set_base_level(*domain);

    // Initialize base level grid
    fclaw2d_domain_iterate_level(*domain, minlevel,
                                 (fclaw2d_patch_callback_t) cb_amrinit,
                                 (void *) NULL);

    // Exchange boundary data at the base level
    bc_level_exchange(*domain,minlevel);

    // Refine as needed.
    for (int level = minlevel; level < maxlevel; level++)
    {
        bool level_refined = false;
        fclaw2d_domain_iterate_level(*domain, level,
                                     (fclaw2d_patch_callback_t) cb_tag_patch,
                                     (void *) &level_refined);
        if (level_refined)
        {
            // Rebuild domain
            fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);

            // This is just for fun; remove when it gets annoying.
            fclaw2d_domain_list_adapted (*domain, new_domain, SC_LP_STATISTICS);

            bool init_flag = true;
            fclaw2d_domain_iterate_adapted(*domain, new_domain,cb_domain_adapt,(void *) &init_flag);

            // Exchange boundary data at the new level
            init_domain_data(new_domain);
            set_domain_data(new_domain, gparms, amropts);
            set_domain_time(new_domain,t);

            // Assume only one block, since we are assuming mthbc
            fclaw2d_block_t *block = &new_domain->blocks[0];
            init_block_data(block);
            set_block_data(block,gparms->m_mthbc);

            int new_level = level+1;
            bc_level_exchange(new_domain,new_level);

            // Deallocate old domain
            amrreset(*domain);
            fclaw2d_domain_destroy(*domain);
            *domain = new_domain;
        }
        else
        {
            // exit loop;  we are done refining
            break;
        }
    }
    cout << "domain (amrinit, end) " << domain << endl;
}


// -----------------------------------------------------------------
// Run - with or without subcycling
// -----------------------------------------------------------------
void amrrun(fclaw2d_domain_t *domain)
{

    cout << "domain (amrrun) " << domain << endl;

    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    cout << "ddata = " << ddata << endl;

    // Write out an initial time file
    int iframe = 0;
    amrout(domain,iframe);

    global_parms *gparms = ddata->parms;
    const amr_options_t *amropts = ddata->amropts;

    // Real final_time = gparms->m_tfinal;
    Real final_time = amropts->tfinal;
    int nout = gparms->m_nout;
    Real initial_dt = gparms->m_initial_dt;

    Real t0 = 0;

    Real dt_outer = (final_time-t0)/Real(nout);
    Real dt_level0 = initial_dt;
    Real t_curr = t0;
    for(int n = 0; n < nout; n++)
    {
        Real tstart = t_curr;
        Real tend = tstart + dt_outer;
        int n_inner = 0;
        while (t_curr < tend)
        {
            subcycle_manager time_stepper;
            time_stepper.define(domain,gparms,ddata->amropts,t_curr);

            // In case we have to reject this step
            save_time_step(domain);
            // check_conservation(domain);

            // Take a stable level 0 time step (use this as the base level time step even if
            // we have no grids on level 0) and reduce it.
            int reduce_factor;
            if (time_stepper.nosubcycle())
            {
                // Take one step of a stable time step for the finest non-emtpy level.
                reduce_factor = time_stepper.maxlevel_factor();
            }
            else
            {
                // Take one step of a stable time step for the coarsest non-empty level.
                reduce_factor = time_stepper.minlevel_factor();
            }
            Real dt_minlevel = dt_level0/reduce_factor;

            // Use the tolerance to make sure we don't take a tiny time step just to
            // hit 'tend'.   We will take a slightly larger time step now (dt_cfl + tol)
            // rather than taking a time step of 'dt_minlevel' now, followed a time step of only
            // 'tol' in the next step.
            // Of course if 'tend - t_curr > dt_minlevel', then dt_minlevel doesn't change.
            Real tol = 1e-2*dt_minlevel;
            bool took_small_step = false;
            if (tend - t_curr - dt_minlevel < tol)
            {
                dt_minlevel = tend - t_curr;  // <= 'dt_minlevel + tol'
                took_small_step = true;
            }

            // This also sets the time step on all finer levels.
            time_stepper.set_dt_minlevel(dt_minlevel);

            Real maxcfl_step = advance_all_levels(domain, &time_stepper);

            printf("Level %d step %5d : dt = %12.3e; maxcfl (step) = %8.3f; Final time = %12.4f\n",
                   time_stepper.minlevel(),n_inner,dt_minlevel,maxcfl_step, t_curr);

            if (maxcfl_step > gparms->m_max_cfl)
            {
                printf("   WARNING : Maximum CFL exceeded; retaking time step\n");
                restore_time_step(domain);

                // Modify dt_level0 from step used.
                dt_level0 = dt_level0*gparms->m_desired_cfl/maxcfl_step;

                // Got back to start of loop, without incrementing step counter or time level
                continue;
            }

            t_curr += dt_minlevel;

            if (took_small_step)
            {
                Real dt0 =  dt_minlevel*reduce_factor;
                printf("   WARNING : Took small time step which was %6.1f%% of desired dt.\n",100.0*dt0/dt_level0);
            }

            // New time step, which should give a cfl close to the desired cfl.
            Real dt_new = dt_level0*gparms->m_desired_cfl/maxcfl_step;
            if (!took_small_step)
            {
                dt_level0 = dt_new;
            }
            else
            {
                // use time step that would have been used had we not taken a small step
            }
            n_inner++;

            int regrid_step = 4;  // Will eventually read this in as a parameter.
            if (n_inner % regrid_step == 0)
            {
                // After some number of time steps, we probably need to regrid...
                // regrid(domain);
            }
        }

        // Output file at every outer loop iteration
        set_domain_time(domain,t_curr);
        iframe++;
        amrout(domain,iframe);
    }
}

void cb_amrout(fclaw2d_domain_t *domain,
               fclaw2d_patch_t *this_patch,
               int this_block_idx,
               int this_patch_idx,
               void *user)
{
    int iframe = *((int *) user);
    int num = this_block_idx*domain->num_blocks + this_patch_idx + 1;
    int matlab_level = this_patch->level + 1;  // Matlab wants levels to start at 1.

    // Patch data is appended to fort.qXXXX
    ClawPatch *cp = get_patch_data(this_patch);

    cp->write_patch_data(iframe, num, matlab_level);
}


void amrout(fclaw2d_domain_t *domain, int iframe)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    cout << "here ....1 " << endl;
    cout << "ddata = " << ddata << endl;
    global_parms *gparms = ddata->parms;
    cout << "here ....1 " << endl;
    Real time = get_domain_time(domain);

    // Get total number of patches
    int ngrids = 0;
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        ngrids += block->num_patches;
    }

    printf("Matlab output Frame %d  at time %12.4f\n\n",iframe,time);

    // Write out header file containing global information for 'iframe'
    write_tfile_(&iframe,&time,&gparms->m_meqn,&ngrids,&gparms->m_maux);

    // This opens file 'fort.qXXXX' for replace (where XXXX = <zero padding><iframe>, e.g. 0001,
    // 0010, 0114), and closes the file.
    new_qfile_(&iframe);

    fclaw2d_domain_iterate_patches(domain, cb_amrout, (void *) &iframe);
}

void amrreset(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *dd;
    dd = (fclaw2d_domain_data_t *) domain->user;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        fclaw2d_block_data_t *bd = (fclaw2d_block_data_t *) block->user;

        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_data_t *pd = (fclaw2d_patch_data_t *) patch->user;
            ClawPatch *cp = pd->cp;

            delete cp;
            FCLAW2D_FREE (pd);
            patch->user = NULL;
        }

        FCLAW2D_FREE (bd);
        block->user = NULL;
    }
    FCLAW2D_FREE (dd);
    domain->user = NULL;
}
