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
using namespace std;

//using std::ifstream;
using std::ios;


void set_domain_data(fclaw2d_domain_t *domain, global_parms *parms)
{
    fclaw2d_domain_data_t *ddata = P4EST_ALLOC (fclaw2d_domain_data_t, 1);
    domain->user = (void *) ddata;

    ddata->parms = parms;
}

global_parms* get_domain_data(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata;
    ddata = (fclaw2d_domain_data_t *) domain->user;

    return ddata->parms;
}

void set_domain_time(fclaw2d_domain_t *domain, Real time)
{
    fclaw2d_domain_data_t *ddata;
    ddata = (fclaw2d_domain_data_t *) domain->user;
    ddata->curr_time = time;
}

Real  get_domain_time(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata;
    ddata = (fclaw2d_domain_data_t *) domain->user;
    return ddata->curr_time;
}


void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp)
{
    fclaw2d_patch_data_t *pdata = P4EST_ALLOC (fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;

    pdata->cp = cp;
}

ClawPatch* get_patch_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata;
    pdata = (fclaw2d_patch_data_t *) patch->user;

    return pdata->cp;
}


// This needs to be defined by p4est
void get_edge_neighbors(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        int iside,
                        int *neighbor_block_idx,
                        int *neighbor_patch_idx,
                        int *relative_refratio)
{
    // Arguments :
    // this_block_idx, this_patch_idx         : Set to the (i,j) indices of the block/patch
    //                                          currently visited by the patch/block loop.
    // iside                                  : the side whose neighbor is requested.  Numbering is :
    //                                             (0,1,2,3)=(left,right,bottom,top).
    // neighbor_block_idx, neighbor_patch_idx : (i,j) index of block and patch corresponding to neighbor at 'iside'
    // relative_refratio                      : Ratio of refinement of neighbor to refinement of
    //                                          'this_patch'.  See below for possible values.
    //

    // Only one patch for now, so there are no neighbors
    *relative_refratio = -1;
}

void get_corner_neighbor(fclaw2d_domain_t *domain,
                         int this_block_idx,
                         int this_patch_idx,
                         int icorner,
                         int *corner_block_idx,
                         int *corner_patch_idx,
                         int *relative_refratio)
{
    // icorner  : corner of 'this_patch' for which neighboring corner is requested. Numbered (0,1,2,3)=(ll,lr,ul,ur)
    //
    // This needs to be defined by p4est
    *relative_refratio = -1;
}

void get_phys_boundary(fclaw2d_domain_t *domain,
                       int this_block_idx,
                       int this_patch_idx,
                       bool *intersects_bc)
{
    // intersects_bc[2*SpaceDim]
    // The entry intersects_bc[iedge]='true' if 'iedge' is at a physical
    // boundary and 'false' otherwise.

    for(int i = 0; i < 2*SpaceDim; i++)
    {
        // Again, only works for runs with a single patch.
        intersects_bc[i] = true;
    }
}


// Called from within fclaw2d_domain_iterate_level
void cb_bc_level_exchange(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    global_parms *gparms = get_domain_data(domain);
    int refratio = gparms->m_refratio;
    ClawPatch *this_cp = get_patch_data(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        int iside = 2*idir + 1;  // Look only at high side for level copy.

        // Output arguments
        int neighbor_block_idx;
        int neighbor_patch_idx[refratio];  // Be prepared to store 1 or more patch indices.
        int ref_flag; // = -1, 0, 1
        get_edge_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iside,
                           &neighbor_block_idx,
                           neighbor_patch_idx,
                           &ref_flag);

        if (ref_flag == 0)
        {
            // Copy data from patch at the same level;
            int num_neighbors = 1; // Same level
            fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
            fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[0]];
            ClawPatch *neighbor_cp = get_patch_data(neighbor_patch);

            // Finally, do exchange between 'this_patch' and 'neighbor patch(es)'.
            this_cp->edge_exchange(idir,iside,num_neighbors,&neighbor_cp);

            // Now check for any physical boundary conditions
            bool intersects_bc[2*SpaceDim];
            Real curr_time = get_domain_time(domain);
            Real dt = 1e20;   // When do we need dt in setting a boundary condition?
            get_phys_boundary(domain,this_block_idx,this_patch_idx,intersects_bc);
            this_cp->set_physbc(intersects_bc,gparms->m_mthbc,curr_time,dt);
        }
    } // loop over directions (idir = 0,1,2)
}


void bc_level_exchange(fclaw2d_domain_t *domain, const int& a_level)
{
    int user = NULL;
    fclaw2d_domain_iterate_level(domain, a_level,
                                  (fclaw2d_patch_callback_t) cb_bc_level_exchange, (void *) user);
}

void cb_bc_average(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    global_parms *gparms = get_domain_data(domain);
    int refratio = gparms->m_refratio;

    ClawPatch *this_cp = get_patch_data(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iside = 2*idir; iside < 2*idir + 1; iside++)
        {
            int neighbor_block_idx;
            int neighbor_patch_idx[refratio];  // Be prepared to store 1 or more patch indices.
            int ref_flag; // = -1, 0, 1
            get_edge_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iside,
                               &neighbor_block_idx,
                               neighbor_patch_idx,
                               &ref_flag);

            if (ref_flag == 1)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'this_patch' by averaging data from finer neighbors
                int num_neighbors = refratio;
                fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                ClawPatch *neighbor_cp[refratio];
                for (int ir = 0; ir < refratio; ir++)
                {
                    fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[ir]];
                    neighbor_cp[ir] = get_patch_data(neighbor_patch);
                }

                // Average finer grid data to coarser 'this_cp' ghost cells.
                // This uses the same routine as 'bc_exchange_level', but now
                // 'num_neighbors = refratio'
                this_cp->edge_exchange(idir,iside,num_neighbors,neighbor_cp);
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
    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    global_parms *gparms = get_domain_data(domain);
    int refratio = gparms->m_refratio;

    ClawPatch *this_cp = get_patch_data(this_patch);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        // Loop over low side and high side
        for (int iside = 2*idir; iside < 2*idir + 1; iside++)
        {
            int neighbor_block_idx;
            int neighbor_patch_idx[refratio];  // Be prepared to store 1 or more patch indices.
            int ref_flag; // = -1, 0, 1
            get_edge_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iside,
                               &neighbor_block_idx,
                               neighbor_patch_idx,
                               &ref_flag);

            if (ref_flag == 1)  // neighbors are at finer level
            {
                // Fill in ghost cells on 'neighbor_patch' by interpolating to finer grid
                int num_neighbors = refratio;
                fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                ClawPatch *neighbor_cp[refratio];
                for (int ir = 0; ir < refratio; ir++)
                {
                    fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[ir]];
                    neighbor_cp[ir] = get_patch_data(neighbor_patch);
                }

                // Average finer grid data to coarser 'this_cp' ghost cells.
                // This uses the same routine as 'bc_exchange_level', but now
                // 'num_neighbors = refratio'
                this_cp->edge_interpolate(idir,iside,num_neighbors,neighbor_cp);
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}


void bc_coarse_exchange(fclaw2d_domain_t *domain, const int& a_level, const int& a_at_time)
{
    // First, average fine grid to coarse grid cells
    int user = NULL;
    int coarser_level = a_level - 1;
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_average,
                                 (void *) user);

    // Interpolate in time, if necessary

    // Interpolate coarse grid to fine.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_interpolate,
                                 (void *) user);
}


void cb_advance_patch(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    fclaw2d_block_t *block = &domain->blocks[this_block_idx];
    fclaw2d_patch_t *patch = &block->patches[this_patch_idx];
    ClawPatch *cp = get_patch_data(patch);
    fclaw2d_level_time_data_t *time_data = (fclaw2d_level_time_data_t *) user;

    Real dt = time_data->dt;
    Real t = time_data->t;

    global_parms *gparms = get_domain_data(domain);
    int refratio = gparms->m_refratio;
    int level = this_patch->level;
    Real maxcfl_grid = cp->step_noqad(t,dt,refratio,level,*gparms);

    //Real maxcfl_grid = cp->step(dt);
    time_data->maxcfl = max(maxcfl_grid,time_data->maxcfl);
}



Real advance_level(fclaw2d_domain_t *domain,
                   const int& a_level,
                   const int& a_from_step,
                   subcycle_manager& a_time_stepper)
{

    // printf("Advancing level %d by dt = %16.4e\n",a_level,a_time_stepper.get_dt(a_level));
    // Check BCs
    Real maxcfl_coarse = 0;
    if (!a_time_stepper.can_advance(a_level,a_from_step))
    {
        if (!a_time_stepper.level_exchange_done(a_level,a_from_step))
        {
            printf("Error (advance_level) : Level exchange at level %d not done at time step %d\n",a_level,a_from_step);
            exit(1);
        }
        if (!a_time_stepper.coarse_exchange_done(a_level,a_from_step))
        {
            int last_coarse_step = a_time_stepper.get_last_step(a_level-1);
            if (a_from_step > last_coarse_step)
            {
                int coarse_level = a_level-1;
                maxcfl_coarse = advance_level(domain,coarse_level,last_coarse_step,a_time_stepper);
            }

            // Level 'a_level' grid will be averaged onto coarser grid ghost cells;  Coarser
            // 'a_level-1' will interpolate to finer grid ghost.
            // This may also require interpolation in time.
            bc_coarse_exchange(domain,a_level,a_from_step);
            int new_coarse_time = a_time_stepper.get_last_step(a_level-1);
            a_time_stepper.set_coarse_exchange(a_level,new_coarse_time);
        }
    }

    fclaw2d_level_time_data_t time_data;

    time_data.maxcfl = maxcfl_coarse;
    time_data.dt = a_time_stepper.get_dt(a_level);
    time_data.t = a_time_stepper.current_time(a_level);

    // Advance this level from 'a_from_time' to 'a_from_time + a_time_stepper.time_step_inc(a_level)'
    fclaw2d_domain_iterate_level(domain, a_level,
                                 (fclaw2d_patch_callback_t) cb_advance_patch,
                                 (void *) &time_data);

    a_time_stepper.increment_time_step(a_level);
    bc_level_exchange(domain,a_level);

    int new_time = a_time_stepper.get_last_step(a_level);
    a_time_stepper.set_level_exchange(a_level,new_time);

    return time_data.maxcfl;  // Maximum from level iteration

}


Real advance_all_levels(fclaw2d_domain_t *domain, const Real& dt)
{

    global_parms *gparms = get_domain_data(domain);
    int maxlevel = gparms->m_maxlevel;
    int refratio = gparms->m_refratio;

    int patches_at_level[maxlevel+1];
    for(int level = 0; level <= maxlevel; level++)
    {
        patches_at_level[level] = num_patches(domain,level);
    }


    // Construct time step manager
    Real t_curr = get_domain_time(domain);
    subcycle_manager time_stepper;
    time_stepper.define(maxlevel,dt,refratio,t_curr,patches_at_level);

    // Time step increment on coarse grid (level 0 grid) is equal
    // to the number of time steps we must take on the fine grid.
    int n_fine_steps = time_stepper.time_step_inc(0);

    // Take 'n_fine_steps' on finest level.  Recursively update coarser levels
    // as needed.
    int finest_level = maxlevel;
    Real maxcfl = 0;
    for(int nf = 0; nf < n_fine_steps; nf++)
    {
        Real cfl_level = advance_level(domain,finest_level,nf,time_stepper);
        maxcfl = max(cfl_level,maxcfl);
    }
    return maxcfl;
}

void amrsetup(fclaw2d_domain_t *domain)
{
    // Check that the minimum level we have is consistent with what
    // was in input file.
    global_parms *gparms = get_domain_data(domain);
    gparms->print_inputParams();

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = new ClawPatch();

            // Set stuff from p4est

            cp->define(patch->xlower,
                       patch->ylower,
                       patch->xupper,
                       patch->yupper,
                       gparms);
            set_patch_data(patch,cp);
        }
    }
}

void cb_amrinit(fclaw2d_domain_t *domain,fclaw2d_patch_t *this_patch,
                int this_block_idx, int this_patch_idx, void *user)
{
    global_parms *gparms = get_domain_data(domain);
    ClawPatch *cp = get_patch_data(this_patch);

    cp->initialize();
    if (gparms->m_maux > 0)
    {
        cp->setAuxArray(gparms->m_maxlevel,gparms->m_refratio,this_patch->level);
    }
}

void amrinit(fclaw2d_domain_t *domain)
{
    Real t = 0;
    set_domain_time(domain,t);

    global_parms *gparms = get_domain_data(domain);
    int minlevel = gparms->m_minlevel;
    int maxlevel = gparms->m_maxlevel;

    // Set problem dependent parameters for Riemann solvers, etc.
    // Values are typically stored in Fortran common blocks, and are not
    // available outside of Fortran.
    setprob_();

    // Don't need to iterate level by level, but saves some coding.
    int user = NULL;
    for(int level = minlevel; level <= maxlevel; level++)
    {
        fclaw2d_domain_iterate_level(domain, level,
                                     (fclaw2d_patch_callback_t) cb_amrinit,
                                     (void *) user);
        bc_level_exchange(domain,level);
    }
}


void amrrun(fclaw2d_domain_t *domain)
{
    // Write out an initial time file
    int iframe = 0;
    amrout(domain,iframe);

    global_parms *gparms = get_domain_data(domain);
    Real final_time = gparms->m_tfinal;
    int nout = gparms->m_nout;
    Real t0 = 0; // Should have this set by user.

    Real dt_outer = (final_time-t0)/Real(nout);

    Real initial_dt = gparms->m_initial_dt;
    Real dt_cfl = initial_dt;
    Real t_curr = t0;
    for(int n = 0; n < nout; n++)
    {
        Real tstart = t_curr;
        Real tend = tstart + dt_outer;
        int n_inner = 0;
        while (t_curr < tend)
        {
            // Use the tolerance to make sure we don't take a tiny time step just to
            // hit 'tend'.   We will take a slightly larger time step now (dt_cfl + tol)
            // rather than taking a time step of 'dt_cfl' now, followed a time step of only
            // 'tol' in the next step.
            // Of course if 'tend - t_curr < dt_cfl', then we take that time step.
            Real dt_inner;
            Real tol = 1e-4*dt_cfl;
            if (tend - t_curr - dt_cfl < tol)
            {
                // In this case, 'dt_curr' is only set temporarily to a smaller value
                dt_inner = tend - t_curr;  // <= 'dt_cfl + tol'
            }
            else
            {
                dt_inner = dt_cfl;
            }

            set_domain_time(domain,t_curr);
            Real maxcfl = advance_all_levels(domain, dt_inner);

            t_curr = t_curr + dt_inner;
            printf("Level 0 step %5d : dt = %12.4e; maxcfl = %8.4f; Final time = %12.4f\n",n_inner,dt_inner,maxcfl, t_curr);

            // New time step, based on max cfl number and desired number.
            dt_cfl = dt_inner*gparms->m_desired_cfl/maxcfl;

            n_inner++;

            // After some number of time steps, we probably need to regrid...
        }

        // Output file at every outer loop iteration
        set_domain_time(domain,t_curr);
        iframe = iframe + 1;
        amrout(domain,iframe);
    }
}


void amrout(fclaw2d_domain_t *domain, int iframe)
{
    global_parms *gparms = get_domain_data(domain);
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
    // 0010, 0114).
    new_qfile_(&iframe);
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = get_patch_data(patch);
            int num = i*domain->num_blocks + j + 1;
            int level = patch->level + 1;

            // Patch data is appended to fort.qXXXX
            cp->write_patch_data(iframe, num, level);
        }
    }
}

void amrreset(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *dd;
    dd = (fclaw2d_domain_data_t *) domain->user;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_data_t *pd = (fclaw2d_patch_data_t *) patch->user;
            ClawPatch *cp = pd->cp;

            delete cp;
            P4EST_FREE (pd);
            patch->user = NULL;
        }
    }
    P4EST_FREE (dd);
    domain->user = NULL;
}
