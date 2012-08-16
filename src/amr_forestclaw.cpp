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


void get_face_neighbors(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        int iside,
                        int *neighbor_block_idx,
                        int neighbor_patch_idx[],
                        int *ref_flag,
                        bool *is_phys_bc)


{
    int rproc[P4EST_HALF];
    int rblockno;
    int rpatchno[P4EST_HALF];
    int rfaceno;

    fclaw2d_face_neighbor_t neighbor_type =
        fclaw2d_patch_face_neighbors(domain,
                                     this_block_idx,
                                     this_patch_idx,
                                     iside,
                                     rproc,
                                     &rblockno,
                                     rpatchno,
                                     &rfaceno);


    // neighbor_type is one of :
    // FCLAW2D_FACE_NEIGHBOR_BOUNDARY,
    // FCLAW2D_FACE_NEIGHBOR_HALFSIZE,
    // FCLAW2D_FACE_NEIGHBOR_SAMESIZE,
    // FCLAW2D_FACE_NEIGHBOR_DOUBLESIZE

    global_parms *gparms = get_domain_data(domain);
    int refratio = gparms->m_refratio; // == P4EST_HALF ??

    *neighbor_block_idx = rblockno;

    if (neighbor_type == FCLAW2D_FACE_NEIGHBOR_BOUNDARY)
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
        if (neighbor_type == FCLAW2D_FACE_NEIGHBOR_HALFSIZE)
        {
            // Neighbors are finer grids
            *ref_flag = 1;
            for(int ir = 0; ir < refratio; ir++)
            {
                neighbor_patch_idx[ir] = rpatchno[ir];
            }
        }
        else if (neighbor_type == FCLAW2D_FACE_NEIGHBOR_SAMESIZE)
        {
            // Neighbor is at the same level
            *ref_flag = 0;
            neighbor_patch_idx[0] = rpatchno[0];
        }
        else if (neighbor_type == FCLAW2D_FACE_NEIGHBOR_DOUBLESIZE)
        {
            // Neighbor is a coarser grid
            *ref_flag = -1;
            neighbor_patch_idx[0] = rpatchno[0];
        }
    }
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


    int bdry[2*SpaceDim]; // 2*SpaceDim == P4EST_FACES
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);

    for(int i = 0; i < 2*SpaceDim; i++)
    {
        intersects_bc[i] = bdry[i] > 0;
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
        bool is_phys_bc;
        get_face_neighbors(domain,
                           this_block_idx,
                           this_patch_idx,
                           iside,
                           &neighbor_block_idx,
                           neighbor_patch_idx,
                           &ref_flag,
                           &is_phys_bc);



        if (ref_flag == 0 && !is_phys_bc)
        {
            // Ghost cells overlap patch at the same level, and so we can do a straight copy.
            // int num_neighbors = 1; // Same level
            fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
            fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[0]];
            ClawPatch *neighbor_cp = get_patch_data(neighbor_patch);
            // Finally, do exchange between 'this_patch' and 'neighbor patch(es)'.
            this_cp->level_face_exchange(idir,iside,&neighbor_cp);
        }
    } // loop over directions (idir = 0,1,2)

    // Now check for any physical boundary conditions on this patch
    bool intersects_bc[2*SpaceDim];
    Real curr_time = get_domain_time(domain);
    Real dt = 1e20;   // When do we need dt in setting a boundary condition?
    get_phys_boundary(domain,this_block_idx,this_patch_idx,intersects_bc);
    this_cp->set_physbc(intersects_bc,gparms->m_mthbc,curr_time,dt);
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
                this_cp->face_average(idir,iside,num_neighbors,neighbor_cp);
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


void bc_exchange_with_coarse(fclaw2d_domain_t *domain, const int& a_level, const int& a_at_time)
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
        printf("Error (advance_level) : Advancing coarser grid (not tested)\n");
        exit(1);
        if (!a_time_stepper->exchanged_with_level(a_level))
        {
            printf("Error (advance_level) : Level exchange at level %d not done at time step %d\n",a_level,a_curr_fine_step);
            exit(1);
        }
        if (!a_time_stepper->exchanged_with_coarser(a_level))
        {
            int last_coarse_step = a_time_stepper->last_step(a_level-1);
            if (a_time_stepper->nosubcycle() || (a_curr_fine_step > last_coarse_step))
            {
                // If we are not subcycling, we this is the only update on the coarser grid that
                // we will do.  We don't however, need this update for getting fine grid
                // boundary conditions.
                // If we are subcycling, 'a_curr_fine_step > last_coarse_step' means that we are trying
                // to advance on the fine grid but can't because the coarse grid hasn't yet
                // advanced to a point beyond the current fine grid time.
                // If we had 'a_curr_fine_step == last_coarse_step', then there would be no need
                // to advance the coarse grid this time around.
                maxcfl = advance_level(domain,a_level-1,last_coarse_step,a_time_stepper);
            }

            // Now we can exchange with the coarser grid.
            // Level 'a_level' grid will be averaged onto coarser grid ghost cells;  Coarser
            // 'a_level-1' will interpolate to finer grid ghost.
            // If we are subcycling, this may also require interpolation in time.  If we are not
            // subcycling, the solution on coarse and fine are updated to the same time level.
            bc_exchange_with_coarse(domain,a_level,a_curr_fine_step);
            a_time_stepper->increment_coarse_exchange_counter(a_level);
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
    // step_inc(minlevel) is the number of steps we must take on the finest level to equal one
    // step on the coarsest non-empty level, i.e. minlevel.
    int minlevel = a_time_stepper->minlevel();
    int n_fine_steps = a_time_stepper->step_inc(minlevel); // equal 1 in the non-subcycled case.

    // Take 'n_fine_steps' on finest level.  Coarser levels will be updated recursively as needed.
    int maxlevel = a_time_stepper->maxlevel();
    Real maxcfl = 0;
    for(int nf = 0; nf < n_fine_steps; nf++)
    {
        Real cfl_step = advance_level(domain,maxlevel,nf,a_time_stepper);
        maxcfl = max(cfl_step,maxcfl);
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

        for(int m = 0; m < 2*SpaceDim; m++)
        {
            block->mthbc[m] = gparms->m_mthbc[m];
        }

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
        printf("\n");
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
            time_stepper.define(domain,t_curr);

            // Take a stable level 0 time step and reduce it.
            int reduce_factor;
            if (time_stepper.nosubcycle())
            {
                // We want the time step to be the stable time step for the finest level.
                reduce_factor = time_stepper.maxlevel_factor();
            }
            else
            {
                // stable time step for the coarsest level that currently has grids.
                reduce_factor = time_stepper.minlevel_factor();
            }
            Real dt_minlevel = dt_level0/reduce_factor;

            // Use the tolerance to make sure we don't take a tiny time step just to
            // hit 'tend'.   We will take a slightly larger time step now (dt_cfl + tol)
            // rather than taking a time step of 'dt_minlevel' now, followed a time step of only
            // 'tol' in the next step.
            // Of course if 'tend - t_curr < dt_minlevel', then dt_minlevel doesn't change.
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
            t_curr += dt_minlevel;

            printf("Level %d step %5d : dt = %12.3e; maxcfl (step) = %8.3f; Final time = %12.4f\n",
                   time_stepper.minlevel(),n_inner,dt_minlevel,maxcfl_step, t_curr);

            if (maxcfl_step > gparms->m_max_cfl)
            {
                printf("   WARNING : Maximum CFL exceeded\n");
            }
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

            // After some number of time steps, we probably need to regrid...
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
    int level = this_patch->level + 1;  // Matlab wants levels to start at 1.

    // Patch data is appended to fort.qXXXX
    ClawPatch *cp = get_patch_data(this_patch);

    cp->write_patch_data(iframe, num, level);
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
