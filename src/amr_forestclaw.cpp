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
#include <cmath>

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

void set_domain_time(fclaw2d_domain_t *domain, double time)
{
    fclaw2d_domain_data_t *ddata;
    ddata = (fclaw2d_domain_data_t *) domain->user;
    ddata->curr_time = time;
}

double get_domain_time(fclaw2d_domain_t *domain)
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
    // This needs to be defined by p4est
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
    // This needs to be defined by p4est
    *relative_refratio = -1;
}

void get_phys_boundary(fclaw2d_domain_t *domain,
                       int this_block_idx,
                       int this_patch_idx,
                       bool *intersects_bc)
{
    for(int i = 0; i < 2*SpaceDim; i++)
    {
        // Again, only works for runs with a single patch.
        intersects_bc[i] = true;
    }
}



void patch_exchange_bc(fclaw2d_domain_t *domain)
{
    // Existing patches should all exchange boundary condition information.
    //
    // Note : Code so far only handles conventional exchanges correctly.
    // More exotic boundaries conditions involving exchanges between, say,
    // two bottom edges (think : two-patch sphere or cubed sphere) or
    // between blocks with different coordinate orientations,
    // (e.g. the mobius strip) are not handled yet. Also, I am trying to
    // write this with 3d in mind, but no guarantees...
    //
    //
    // Steps :
    // Step 1 : Average all fine grids to coarse grid ghost cells;  Exchange values between grids
    //          that are at the same level.  Set corner values where no interpolation is required.
    // Step 2 : Set physical boundary conditions.
    // Step 3 : Interpolate coarse grid data to fine grid ghost cell values.
    //

    global_parms *gparms = get_domain_data(domain);
    int refratio = gparms->m_refratio;



    // Step 1 :
    // Exchange boundary data with neighboring patches that are at 'this_level'
    // or a finer level.  The current patch should inititiate exchange with
    // neighbors.
    //
    // To avoid duplication of work, if the neighbor level is equal
    // to 'this_level', then the current patch will only initiate the exchange
    // with a high side patch. The low side edge will be filled by with
    // exchanges initiated by the low side patch when it is visited in the patch loop.
    //
    // Part of this step also involves exchanging corner ghost cell
    // information at this level or finer levels.
    //
    // Interpolation from coarse to fine grids, as well as physical BCs
    // will be handled in a later step.
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *this_patch = &block->patches[j];
            ClawPatch *this_cp = get_patch_data(this_patch);
            // int this_level = this_patch->level;
            int this_block_idx = i;
            int this_patch_idx = j;

            for (int idir = 0; idir < SpaceDim; idir++)
            {
                int neighbor_block_idx;
                int neighbor_patch_idx[refratio];  // Be prepared to store 1 or more patch indices.
                int relative_refratio;

                // 'relative_refratio' is the refinement ratio of neighbor patch(es) to 'this_patch'
                //
                // -1           Neighbor is not valid for this step (either because it is a
                //                   coarser patch, the boundary is a physical boundary, or the
                //                   neighbor is a low side patch at 'this_level')
                // 1            High-side patch is at 'this_level'
                // 'refratio'   Patch has 'refratio' finer grids.
                //

                // Loop over low side then high side
                for (int iside = 2*idir; iside < 2*idir + 1; iside++)
                {
                    get_edge_neighbors(domain,
                                       this_block_idx,
                                       this_patch_idx,
                                       iside,
                                       &neighbor_block_idx, // Am I passing a pointer?
                                       neighbor_patch_idx,
                                       &relative_refratio);

                    if (relative_refratio > 0)  // check for valid neighbor patches
                    {
                        fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                        fclaw2d_patch_t *neighbor_patch[relative_refratio];
                        ClawPatch *neighbor_cp[relative_refratio];
                        for (int ir = 0; ir < relative_refratio; ir++)
                        {
                            neighbor_patch[ir]  = &neighbor_block->patches[neighbor_patch_idx[ir]];
                            neighbor_cp[ir] = get_patch_data(neighbor_patch[ir]);
                        }

                        this_cp->edge_exchange_step1(iside,relative_refratio,neighbor_cp);
                    } // relative_refratio

                    // Exchange corner ghost cells between patches at the same level
                    // Corners ordered (0=ll, 1=lr, 2=ul, 3=ur)

                    int corner_block_idx;
                    int corner_patch_idx;
                    int num_corners = pow(SpaceDim,2);
                    for (int icorner = 0; icorner < num_corners; icorner++)
                    {
                        get_corner_neighbor(domain,
                                            this_block_idx,
                                            this_patch_idx,
                                            icorner,
                                            &corner_block_idx,
                                            &corner_patch_idx,
                                            &relative_refratio);
                        if (relative_refratio >= 0)
                        {
                            fclaw2d_block_t *corner_block = &domain->blocks[corner_block_idx];
                            fclaw2d_patch_t *corner_patch = &corner_block->patches[corner_patch_idx];
                            ClawPatch *corner_cp = get_patch_data(corner_patch);
                            this_cp->corner_exchange_step1(icorner,relative_refratio,corner_cp);
                        }
                    } // loop over corners
                } // loop over sides (lo --> hi)
            } // loop over directions (idir = 0,1,2)
        } // loop over patches on block
    } // loop over blocks in domain


    // Step 2 :
    // Set physical boundary conditions.
    Real curr_time = get_domain_time(domain);
    Real dt = 1e20;   // I am not sure how this is used in setting boundary conditions
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *this_patch = &block->patches[j];
            ClawPatch *this_cp = get_patch_data(this_patch);
            int this_block_idx = i;
            int this_patch_idx = j;

            bool intersects_bc[2*SpaceDim];
            get_phys_boundary(domain,this_block_idx,this_patch_idx,intersects_bc);
            this_cp->set_physbc_step2(intersects_bc,gparms->m_mthbc,curr_time,dt);
        }
    }

    // I should probably now handle the few cases where  neighboring patches share a common
    // physical boundary.  In this case, corner ghost cells from each patch will overlap ghost
    // cells from the neighboring patch.   These corner cells should be properly filled in....
    // ....


    // Step3 : Interpolate coarse grid data to fine grid ghost cells.
    // Set physical boundary conditions.
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *this_patch = &block->patches[j];
            ClawPatch *this_cp = get_patch_data(this_patch);
            int this_block_idx = i;
            int this_patch_idx = j;
            for (int idir = 0; idir < SpaceDim; idir++)
            {
                int neighbor_block_idx;
                int neighbor_patch_idx[refratio];  // Be prepared to store 1 or more patch indices.
                int relative_refratio;

                // Loop over low side then high side
                for (int iside = 2*idir; iside < 2*idir + 1; iside++)
                {
                    get_edge_neighbors(domain,
                                       this_block_idx,
                                       this_patch_idx,
                                       iside,
                                       &neighbor_block_idx, // Am I passing a pointer?
                                       neighbor_patch_idx,
                                       &relative_refratio);

                    if (relative_refratio == refratio)  // check that we have fine grid neighbors
                    {
                        fclaw2d_block_t *neighbor_block = &domain->blocks[neighbor_block_idx];
                        fclaw2d_patch_t *neighbor_patch[relative_refratio];
                        ClawPatch *neighbor_cp[relative_refratio];
                        for (int ir = 0; ir < relative_refratio; ir++)
                        {
                            neighbor_patch[ir]  = &neighbor_block->patches[neighbor_patch_idx[ir]];
                            neighbor_cp[ir] = get_patch_data(neighbor_patch[ir]);
                        }

                        this_cp->edge_exchange_step3(iside,refratio,neighbor_cp);
                    } // relative_refratio
                }
            }
        }
    }
}


void amrsetup(fclaw2d_domain_t *domain)
{
    global_parms *gparms = new global_parms();
    cout << "Global parameters " << endl;
    gparms->get_inputParams();
    gparms->print_inputParams();

    set_domain_data(domain,gparms);

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

void amrinit(fclaw2d_domain_t *domain)
{
    double t = 0;
    set_domain_time(domain,t);

    global_parms *gparms = get_domain_data(domain);

    // Set problem dependent parameters for Riemann solvers, etc.
    // Values are typically stored in Fortran common blocks, and are not
    // available outside of Fortran.
    setprob_();

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = get_patch_data(patch);

            cp->initialize();
            cp->setAuxArray(gparms->m_maxlevel,gparms->m_refratio,patch->level);
        }
    }
}


void amrrun(fclaw2d_domain_t *domain)
{
    // Write out an initial time file
    int iframe = 0;
    amrout(domain,iframe);

    // global_parms *gparms = get_domain_data(domain);
    // int refratio = gparms->m_refratio;

    // Do fake timestepping for now
    for(int iframe = 1; iframe < 5; iframe++)
    {
        // First exchange boundary conditions
        patch_exchange_bc(domain);
        for(int i = 0; i < domain->num_blocks; i++)
        {
            fclaw2d_block_t *block = &domain->blocks[i];
            for(int j = 0; j < block->num_patches; j++)
            {
                // Fake update
                cout << "Updating solution on patch number " << j << endl;
            }
        }
        amrout(domain,iframe);
    }
}


void amrout(fclaw2d_domain_t *domain, int iframe)
{
    global_parms *gparms = get_domain_data(domain);
    double time = get_domain_time(domain);

    // Get total number of patches
    int ngrids = 0;
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        ngrids += block->num_patches;
    }

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
