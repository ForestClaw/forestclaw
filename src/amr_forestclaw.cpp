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

class ClawPatch;
class global_parms;

void explicit_step_fixed_output(fclaw2d_domain_t *domain);
void explicit_step(fclaw2d_domain_t *domain, int nstep, int nplot);

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

// -----------------------------------------------------------------
// Initialize data
// -----------------------------------------------------------------
static void init_domain_data(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = FCLAW2D_ALLOC_ZERO (fclaw2d_domain_data_t, 1);
    domain->user = (void *) ddata;
}

static void init_block_data(fclaw2d_block_t *block)
{
    fclaw2d_block_data_t *bdata = FCLAW2D_ALLOC_ZERO (fclaw2d_block_data_t, 1);
    block->user = (void *) bdata;
}


static void init_patch_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;
}

// -----------------------------------------------------------------
// Return pointer to user data
// -----------------------------------------------------------------
fclaw2d_domain_data_t *get_domain_data(fclaw2d_domain_t *domain)
{
    return (fclaw2d_domain_data_t *) domain->user;
}


fclaw2d_block_data_t *get_block_data(fclaw2d_block_t *block)
{
    return (fclaw2d_block_data_t *) block->user;
}


fclaw2d_patch_data_t *get_patch_data(fclaw2d_patch_t *patch)
{
    return (fclaw2d_patch_data_t *) patch->user;
}


// -----------------------------------------------------------------
// Set user data with user defined variables, etc.
// -----------------------------------------------------------------
void set_domain_data(fclaw2d_domain_t *domain, global_parms *gparms, const amr_options_t *amropts)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    ddata->parms = gparms;
    ddata->amropts = amropts;
}

void set_block_data(fclaw2d_block_t *block, int *mthbc)
{
    fclaw2d_block_data_t *bdata = get_block_data(block);
    for(int i = 0; i < 4; i++)
    {
        bdata->mthbc[i] = mthbc[i];
    }
}

void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(patch);
    pdata->cp = cp;
}


// -----------------------------------------------------------------
// Some lazy helper functions that really do make things easier..
// -----------------------------------------------------------------
void allocate_user_data(fclaw2d_domain_t *domain)
{
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    init_domain_data(domain);

    for (int i = 0; i < domain->num_blocks; i++)
    {
        block = &domain->blocks[i];
        init_block_data(block);
        for (int j = 0; j < block->num_patches; j++)
        {
            patch = &block->patches[j];
            init_patch_data(patch);
        }
    }
}


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
ClawPatch* get_clawpatch(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t *) patch->user;

    return pdata->cp;
}
/* end of helper functions */

const int get_refratio(fclaw2d_domain_t *domain)
{
    global_parms *gparms = get_domain_parms(domain);
    return gparms->m_refratio;
}

const int get_corners_per_patch(fclaw2d_domain_t *domain)
{
    // Number of patch corners, not the number of corners in the domain!
    return fclaw2d_domain_num_corners(domain);
}

const int get_faces_per_patch(fclaw2d_domain_t *domain)
{
    // Number of faces per patch, not the total number of faces in the domain!
    return fclaw2d_domain_num_faces(domain);
}


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
    ClawPatch *this_cp = get_clawpatch(this_patch);
    Real *sum = (Real*) user;

    *sum += this_cp->compute_sum();
}


void check_conservation(fclaw2d_domain_t *domain)
{
    Real sum = 0;
    fclaw2d_domain_iterate_patches(domain,cb_check_conservation,(void *) &sum);

    printf("Total sum = %24.16f\n",sum);
}

// Dump current patch
void cb_dump_patch(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patch = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patch + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump();
    }
}

void dump_patch(fclaw2d_domain_t *domain, int dump_patch)
{
    printf("Dumping patch (current) %d\n",dump_patch);
    fclaw2d_domain_iterate_patches(domain,(fclaw2d_patch_callback_t) cb_dump_patch,
                                   &dump_patch);
}

// Dump last patch
void cb_dump_last_patch(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patch = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patch + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump_last();
    }
}

void dump_last_patch(fclaw2d_domain_t *domain, int dump_patch)
{
    printf("Dumping patch (last) %d\n",dump_patch);
    fclaw2d_domain_iterate_patches(domain,(fclaw2d_patch_callback_t) cb_dump_last_patch,
                                   &dump_patch);
}
void cb_dump_time_interp_patch(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patch = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patch + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump_time_interp();
    }
}

void dump_time_interp_patch(fclaw2d_domain_t *domain, int dump_patch)
{
    printf("Dumping patch (time_interp) %d\n",dump_patch);
    fclaw2d_domain_iterate_patches(domain,
                                   (fclaw2d_patch_callback_t) cb_dump_time_interp_patch,
                                   &dump_patch);
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
                        int **ref_flag_ptr)
{
    const int p4est_refineFactor = fclaw2d_domain_num_face_corners (domain);
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
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
        // Set the pointer to NULL rather than come up with some bogus value for ref_flag.
        *ref_flag_ptr = NULL;
     }
    else
    {
        if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            // Neighbors are finer grids
            **ref_flag_ptr = 1;
            for(int ir = 0; ir < p4est_refineFactor; ir++)
            {
                neighbor_patch_idx[ir] = rpatchno[ir];
            }
        }
        else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            // Neighbor is at the same level
            **ref_flag_ptr = 0;
            neighbor_patch_idx[0] = rpatchno[0];
        }
        else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
        {
            // Neighbor is a coarser grid
            **ref_flag_ptr = -1;
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
                         int **ref_flag_ptr)
{
    fclaw2d_patch_relation_t neighbor_type;
    int rproc;
    int has_corner_neighbor =
      fclaw2d_patch_corner_neighbors (domain, this_block_idx, this_patch_idx,
        icorner, &rproc, corner_block_idx, corner_patch_idx, &neighbor_type);

    if (!has_corner_neighbor)
    {
        // printf("get_corner_neighbors : Patch %d at corner %d does not have any "
        //        "corner neighbors\n",this_patch_idx,icorner);
        *ref_flag_ptr = NULL;
    }
    else if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
    {
      **ref_flag_ptr = 1;
    }
    else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
    {
      **ref_flag_ptr = 0;
    }
    else
    {
      **ref_flag_ptr = -1;
    }
}


void get_block_boundary(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        bool *intersects_block)
{
    const int p4est_refineFactor = fclaw2d_domain_num_face_corners (domain);
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;
    const int numfaces = get_faces_per_patch(domain);

    for (int iside = 0; iside < numfaces; iside++)
    {
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
        if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
        {
            // 'iside' is a physical boundary
            intersects_block[iside] = true;
        }
        else
        {
            // We have a neighbor patch on block 'rblockno'.
            intersects_block[iside] = this_block_idx != rblockno;
        }
    }
}

void get_phys_boundary(fclaw2d_domain_t *domain,
                       int this_block_idx,
                       int this_patch_idx,
                       bool *intersects_bc)
{
    const int numfaces = get_faces_per_patch(domain);
    int bdry[numfaces];
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);
    for(int i = 0; i < numfaces; i++)
    {
        // Physical boundary conditions
        intersects_bc[i] = bdry[i] == 1;
    }
}


// -----------------------------------------------------------------
// Exchange corner and face information at same level
// -----------------------------------------------------------------
void cb_bc_level_face_exchange(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               void *user)
{
    const int p4est_refineFactor = fclaw2d_domain_num_face_corners(domain);
    // const int numfaces = get_faces_per_patch(domain);
    ClawPatch *this_cp = get_clawpatch(this_patch);

    int numfaces = get_faces_per_patch(domain);
    for (int iface = 0; iface < numfaces; iface++)
    {
        // Output arguments
        int neighbor_block_idx;
        int neighbor_patch_idx[p4est_refineFactor];  // Be prepared to store 1 or more patch indices.
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
                    if (this_patch_idx == 3 && this_patch->level == 3)
                    {
                        // cout << "this_neighbor_patch = " << this_neighbor_idx << endl;
                        // dump_patch(domain,3);
                    }

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
        }
    } // loop over faces
}

void cb_level_corner_exchange(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              void *user)
{

    const int numfaces = get_faces_per_patch(domain);
    bool intersects_bc[numfaces];
    bool intersects_block[numfaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    // Number of patch corners, not the number of corners in the domain!
    const int numcorners = get_corners_per_patch(domain);

    for (int icorner = 0; icorner < numcorners; icorner++)
    {
        // p4est has tons of lookup table like this, can be exposed similarly
        int corner_faces[SpaceDim];
        fclaw2d_domain_corner_faces(domain, icorner, corner_faces);

        // Both faces are at a physical boundary
        bool is_phys_corner =
                intersects_bc[corner_faces[0]] && intersects_bc[corner_faces[1]];

        // Corner lies in interior of physical boundary edge.
        bool corner_on_phys_face = !is_phys_corner &&
                (intersects_bc[corner_faces[0]] || intersects_bc[corner_faces[1]]);

        // Either a corner at a block boundary (but not a physical boundary),
        // or internal to a block
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
            // Both faces are at a physical boundary
            bool is_block_corner =
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


// this routine can probably be combined with cb_level_corner_exchange
void cb_corner_average(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    const int numfaces = get_faces_per_patch(domain);
    const int refratio = get_refratio(domain);
    bool intersects_bc[numfaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    // Number of patch corners, not the number of corners in the domain!
    const int numcorners = get_corners_per_patch(domain);

    for (int icorner = 0; icorner < numcorners; icorner++)
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

                fclaw2d_subcycle_info_t *step_info = (fclaw2d_subcycle_info_t*) user;
                bool time_interp = step_info->do_time_interp;

                // 'this' is the coarser level; 'corner_cp' is the finer level
                this_cp->average_corner_ghost(icorner,refratio,corner_cp,time_interp);
            }
        }
    }
}

// this routine can probably be combined with cb_level_corner_exchange
void cb_corner_interpolate(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           void *user)
{
    const int numfaces = get_faces_per_patch(domain);
    const int refratio = get_refratio(domain);
    bool intersects_bc[numfaces];

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bc);

    // Number of patch corners, not the number of corners in the domain!
    const int numcorners = get_corners_per_patch(domain);

    for (int icorner = 0; icorner < numcorners; icorner++)
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

                fclaw2d_subcycle_info_t *step_info = (fclaw2d_subcycle_info_t*) user;
                bool time_interp = step_info->do_time_interp;

                // 'this' is the coarser level; 'corner_cp' is the finer level
                this_cp->interpolate_corner_ghost(icorner,refratio,corner_cp,time_interp);
            }
        }
    }
}

void cb_set_phys_bc(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    int numfaces = get_faces_per_patch(domain);
    bool intersects_bc[numfaces];
    Real curr_time = *((Real*) user);
    Real dt = 1e20;   // When do we need dt in setting a boundary condition?
    get_phys_boundary(domain,this_block_idx,this_patch_idx,intersects_bc);

    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    fclaw2d_block_data_t *bdata = get_block_data (this_block);

    ClawPatch *this_cp = get_clawpatch(this_patch);
    this_cp->set_phys_face_ghost(intersects_bc,bdata->mthbc,curr_time,dt);
}

void bc_set_phys(fclaw2d_domain_t *domain, int a_level, Real a_level_time)
{
    fclaw2d_domain_iterate_level(domain, a_level,
                                 (fclaw2d_patch_callback_t) cb_set_phys_bc,
                                 (void *) &a_level_time);
}

void bc_level_exchange(fclaw2d_domain_t *domain, int a_level)
{
    fclaw2d_subcycle_info step_info;
    step_info.level_time = get_domain_time(domain);

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
    const int p4est_refineFactor = fclaw2d_domain_num_face_corners(domain);

    // We may need to average onto a time interpolated grid, not the actual solution.
    fclaw2d_subcycle_info *step_info = (fclaw2d_subcycle_info*) user;
    bool do_time_interp = step_info->do_time_interp;

    // Fill in ghost cells at level 'a_level' by averaging from level 'a_level + 1'
    const int refratio = get_refratio(domain);

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
                    fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[ir]];
                    fine_neighbor_cp[ir] = get_clawpatch(neighbor_patch);
                }
                bool block_boundary = this_block_idx != neighbor_block_idx;
                // Fill in ghost cells on 'this_cp' by averaging from 'fine_neighbor_cp'
                this_cp->average_face_ghost(idir,iface,p4est_refineFactor,refratio,
                                            fine_neighbor_cp,do_time_interp,block_boundary);
            }
        } // loop sides (iside = 0,1,2,3)
    } // loop over directions (idir = 0,1,2)
}


// Iterator over patches looking for finer neighbors
void cb_bc_interpolate(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    const int p4est_refineFactor = fclaw2d_domain_num_face_corners(domain);
    const int refratio = get_refratio(domain);

    fclaw2d_subcycle_info_t *step_info = (fclaw2d_subcycle_info_t*) user;

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
                    fclaw2d_patch_t *neighbor_patch = &neighbor_block->patches[neighbor_patch_idx[ir]];
                    fine_neighbor_cp[ir] = get_clawpatch(neighbor_patch);
                }

                // Fill in fine grid ghost on 'fine_neighbor_cp' by interpolating from 'this_cp',
                // doing time interpolation if necessary
                bool block_boundary = this_block_idx != neighbor_block_idx;
                this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                refratio,fine_neighbor_cp,step_info->do_time_interp,
                                                block_boundary);
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

    ClawPatch *cp = get_clawpatch(this_patch);
    cp->time_interpolate(step_info->fine_step, step_info->coarse_step, step_info->refratio);
}

void bc_exchange_with_coarse_time_interp(fclaw2d_domain_t *domain, const int& a_level,
                                         const int& a_coarse_step, const int& a_fine_step,
                                         const int& a_refratio)
{
    // First, average fine grid to coarse grid cells
    fclaw2d_subcycle_info_t step_info;
    step_info.coarse_step = a_coarse_step;
    step_info.fine_step = a_fine_step;
    step_info.refratio = a_refratio;
    step_info.do_time_interp = true;

    Real level_time = get_domain_time(domain);


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


    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_corner_average,
                                 (void *) &step_info);

    // This is needed so that interpolation below works near boundary.
    bc_set_phys(domain,coarser_level,level_time);

    // Interpolate coarse grid to fine.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_interpolate,
                                 (void *) &step_info);


    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_corner_interpolate,
                                 (void *) &step_info);

}

void bc_exchange_with_coarse(fclaw2d_domain_t *domain, const int& a_level)
{
    // Simple exchange - no time interpolation needed
    fclaw2d_subcycle_info_t step_info;
    Real level_time = get_domain_time(domain);
    step_info.level_time = level_time;
    step_info.do_time_interp = false;

    // Iterate over coarser level and average from finer neighbors to coarse.
    int coarser_level = a_level - 1;
    fclaw2d_domain_iterate_level(domain, coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_average,
                                 (void *) &step_info);

    // Average fine grid corners to the coarse grid ghost cells
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_corner_average,
                                 (void *) &step_info);

    bc_set_phys(domain,coarser_level,level_time);

    // Interpolate coarse grid to fine.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_bc_interpolate,
                                 (void *) &step_info);

    // Interpolate coarse grid to fine grid ghost cells.
    fclaw2d_domain_iterate_level(domain,coarser_level,
                                 (fclaw2d_patch_callback_t) cb_corner_interpolate,
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
    ClawPatch *this_cp = get_clawpatch(this_patch);

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
    ClawPatch *this_cp = get_clawpatch(this_patch);

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

    ClawPatch *cp = get_clawpatch(this_patch);
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
    bool verbose = false;
    Real t_level = a_time_stepper->current_time(a_level);

    Real maxcfl = 0;
    if (verbose)
    {
        cout << endl;
        cout << "Advancing level " << a_level << " from step " <<
            a_curr_fine_step << endl;
    }

    if (!a_time_stepper->can_advance(a_level,a_curr_fine_step))
    {
        if (!a_time_stepper->level_exchange_done(a_level))
        {
            // Level exchange should have been done right after solution update.
            printf("Error (advance_level) : Level exchange at level %d not done at time step %d\n",
                   a_level,a_curr_fine_step);
            exit(1);
        }
        if (!a_time_stepper->exchanged_with_coarser(a_level))
        {
            int last_coarse_step = a_time_stepper->last_step(a_level-1);
            if (verbose)
            {
                cout << " --> Exchange between fine level " << a_level <<
                    " and coarse level " << a_level-1 << " not done at step " << a_curr_fine_step << endl;
            }
            if (a_curr_fine_step == last_coarse_step)
            {
                // Levels are time synchronized and we can do a simple coarse/fine
                // exchange without time interpolation or advancing the coarser level
                if (verbose)
                {
                    cout << " ----> Coarse and fine level are time synchronized; doing exchange"
                        " without time interpolation" << endl;
                }
                bc_exchange_with_coarse(domain,a_level);
                bc_set_phys(domain,a_level,t_level);
                a_time_stepper->increment_coarse_exchange_counter(a_level);
                a_time_stepper->increment_fine_exchange_counter(a_level-1);

                if (a_time_stepper->nosubcycle() && !a_time_stepper->is_coarsest(a_level))
                {
                    // non-subcycled case : this advance is a bit gratuitous, because we
                    // don't need it to advance the fine grid;  rather, we put this advance
                    // here as a way to get advances on the coarser levels.
                    if (verbose)
                    {
                        cout << " ----> Non-subcycled case " << endl;
                        cout << " ----> Making recursive call to advance_level for level " << a_level-1 << endl;
                    }
                    maxcfl = advance_level(domain,a_level-1,last_coarse_step,a_time_stepper);
                    if (verbose)
                    {
                        cout << " ----> Returning from recursive call at level " << a_level << endl;
                    }
                }
            }
            else
            {
                if (verbose)
                {
                    cout << " --> Coarse and fine level are not time synchronized; doing exchange "
                        "with time interpolation" << endl;
                }
                if ((a_curr_fine_step > last_coarse_step))
                {
                    // subcycled case : a_curr_fine_step will only be greater than
                    // last_coarse_step if we haven't yet advanced the coarse grid to a time level
                    // beyond the current fine grid level.
                    if (verbose)
                    {
                        cout << " ----> Subcycled case " << endl;
                        cout << " ----> Making recursive call to advance_level for level " << a_level-1 << endl;
                    }
                    maxcfl = advance_level(domain,a_level-1,last_coarse_step,a_time_stepper);
                    if (verbose)
                    {
                        cout << " ----> Returning from recursive call at level " << a_level << endl;
                    }
                }
                if (!a_time_stepper->nosubcycle())
                {
                    if (verbose)
                    {
                        cout << " --> Doing time interpolatation from coarse grid at level " << a_level-1 << endl;
                    }
                    int refratio = get_refratio(domain);

                    // (1) a_curr_fine_step > last_coarse_step : we just advanced the coarse grid
                    // and can now apply time interpolated boundary conditions, or
                    //
                    // (2) a_curr_fine_step < last_coarse_step : we advanced the coarse
                    // grid in a previous step but we still have to do time interpolation (this
                    // can only happen if refratio > 2)

                    bc_exchange_with_coarse_time_interp(domain,a_level,last_coarse_step,
                                                        a_curr_fine_step,refratio);
                    bc_set_phys(domain,a_level,t_level);
                    a_time_stepper->increment_coarse_exchange_counter(a_level);

                    // Don't increment the fine_exchange_counter, since in the time interpolated case,
                    // there is no coarse time data at time step a_curr_fine_step.
                    // a_time_stepper->increment_fine_exchange_counter(a_level-1);
                }
            }
        }
    }
    if (verbose)
    {
        cout << "Taking step on level " << a_level << endl;
    }

    fclaw2d_level_time_data_t time_data;

    time_data.maxcfl = maxcfl;
    time_data.dt = a_time_stepper->dt(a_level);
    time_data.t = t_level;

    // Advance this level from 'a_curr_fine_step' to 'a_curr_fine_step +
    // a_time_stepper.step_inc(a_level)'
    fclaw2d_domain_iterate_level(domain, a_level,
                                 (fclaw2d_patch_callback_t) cb_advance_patch,
                                 (void *) &time_data);
    a_time_stepper->increment_step_counter(a_level);
    a_time_stepper->increment_time(a_level);

    bc_level_exchange(domain,a_level);

    bc_set_phys(domain,a_level,t_level);
    a_time_stepper->increment_level_exchange_counter(a_level);

    if (verbose)
    {
        cout << "Advance on level " << a_level << " done" << endl << endl;
    }

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

    bool init_flag = *((bool *) user);
    global_parms *gparms = get_domain_parms(domain);
    int maxlevel = gparms->m_maxlevel;
    int minlevel = gparms->m_minlevel;
    int refratio = gparms->m_refratio;
    bool patch_refined = false;
    bool patch_coarsened = false;

    ClawPatch *cp = get_clawpatch(this_patch);

    int level = this_patch->level;
    if (level < maxlevel)
    {
        patch_refined = cp->tag_for_refinement(init_flag);

        if (patch_refined)
        {
            fclaw2d_patch_mark_refine(domain, this_block_idx, this_patch_idx);
        }
    }

    // If a patch has been tagged for refinement, then we shouldn't coarsen it.
    if (level > minlevel && !init_flag && !patch_refined)
    {
        patch_coarsened = cp->tag_for_coarsening(refratio);
        if (patch_coarsened)
        {
            fclaw2d_patch_mark_coarsen(domain, this_block_idx, this_patch_idx);
        }
    }
    if (patch_refined && patch_coarsened)
    {
        printf("Patch tagged for both refinement and coarsening\n");
        exit(1);
    }
}


void cb_init_base(fclaw2d_domain_t *domain,
                  fclaw2d_patch_t *this_patch,
                  int this_block_idx,
                  int this_patch_idx,
                  void *user)
{
    global_parms *gparms = get_domain_parms(domain);
    ClawPatch *cp = new ClawPatch();

    cp->define(this_patch->xlower,
               this_patch->ylower,
               this_patch->xupper,
               this_patch->yupper,
               this_block_idx,
               gparms);

    int level  = this_patch->level;
    int refratio = gparms->m_refratio;
    int maxlevel = gparms->m_maxlevel;

    cp->setup_patch(level, maxlevel, refratio);
    set_patch_data(this_patch,cp);
}


void amr_set_base_level(fclaw2d_domain_t *domain, const int& level)
{
    global_parms *gparms = get_domain_parms(domain);

    fclaw2d_domain_iterate_level(domain, level,
                                 (fclaw2d_patch_callback_t) cb_init_base,
                                 (void *) gparms);
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
    global_parms *gparms = get_domain_parms(old_domain);

    const int num_siblings = fclaw2d_domain_num_corners (old_domain);
    bool init_grid = *(bool *) user;

    int refratio = gparms->m_refratio;
    int maxlevel = gparms->m_maxlevel;

    if (newsize == FCLAW2D_PATCH_SAMESIZE)
    {
        // Grid was not coarsened or refined.
        ClawPatch *cp_old = get_clawpatch(&old_patch[0]);
        ClawPatch *cp_new = new ClawPatch();
        cp_new->define(old_patch[0].xlower,
                       old_patch[0].ylower,
                       old_patch[0].xupper,
                       old_patch[0].yupper,
                       blockno,
                       gparms);

        int level = old_patch[0].level;
        cp_new->copyFrom(cp_old);  // Copy grid data and aux data
        cp_new->setup_patch(level,maxlevel, refratio);
        set_patch_data(&new_patch[0],cp_new);
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        // New grids are FINER grids
        ClawPatch *cp_old = get_clawpatch(&old_patch[0]);

        for (int igrid = 0; igrid < num_siblings; igrid++)
        {
            ClawPatch *cp_new = new ClawPatch();

            cp_new->define(new_patch[igrid].xlower,
                           new_patch[igrid].ylower,
                           new_patch[igrid].xupper,
                           new_patch[igrid].yupper,
                           blockno,
                           gparms);

            int level = new_patch[igrid].level;
            cp_new->setup_patch(level, maxlevel, refratio);
            if (init_grid)
            {
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
                       blockno,
                       gparms);

        int level = new_patch[0].level;
        cp_new->setup_patch(level, maxlevel, refratio);
        for(int igrid = 0; igrid < num_siblings; igrid++)
        {
            // Get one of the older finer grids and average to new coarse grid.
            // Assume that we will never coarsen when we are initializing the grids.
            ClawPatch *cp_old = get_clawpatch(&old_patch[igrid]);
            cp_new->coarsen_from_fine_patch(cp_old, igrid, p4est_refineFactor,refratio);
        }
        set_patch_data(&new_patch[0],cp_new);
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
void cb_amrinit(fclaw2d_domain_t *domain,
                fclaw2d_patch_t *this_patch,
                int this_block_idx,
                int this_patch_idx,
                void *user)
{
    ClawPatch *cp = get_clawpatch(this_patch);

    cp->initialize();
}

// Initialize a base level of grids
void amrinit(fclaw2d_domain_t **domain,
             global_parms *gparms,
             const amr_options_t * amropts)
{
    Real t = 0;

    // gparms->print_inputParams();

    allocate_user_data(*domain);

    set_domain_data(*domain, gparms, amropts);
    set_domain_time(*domain,t);

    int minlevel = gparms->m_minlevel;
    int maxlevel = gparms->m_maxlevel;

    // Set problem dependent parameters for Riemann solvers, etc.
    // Values are typically stored in Fortran common blocks, and are not
    // available outside of Fortran.
    set_problem_parameters();

    // Set up storage for base level grids so we can initialize them
    // Allocates per-block and per-patch user data

    // This function is redundant, and should be made more general.
    cout << "Setting base level " << endl;
    amr_set_base_level(*domain,minlevel);

    cout << "Done with amr_set_base_level " << endl;


    // Initialize base level grid - combine with 'amr_set_base_level' above?
    fclaw2d_domain_iterate_level(*domain, minlevel,
                                 (fclaw2d_patch_callback_t) cb_amrinit,
                                 (void *) NULL);

    cout << "Done with domain adaptation " << endl;

    int num = (*domain)->num_blocks;
    for (int i = 0; i < num; i++)
    {
        fclaw2d_block_t *block = &(*domain)->blocks[i];
        set_block_data(block,gparms->m_mthbc);
    }

    bc_set_phys(*domain,minlevel,t);

    // Refine as needed.

    bool init_flag = true;
    for (int level = minlevel; level < maxlevel; level++)
    {
        cout << "amrinit : level = " << level << endl << endl;

        // TODO: don't use level_refined since it is not agreed upon in parallel
        // the fclaw2d_domain_adapt and _partition calls work fine in parallel

        fclaw2d_domain_iterate_level(*domain, level,
                                     (fclaw2d_patch_callback_t) cb_tag_patch,
                                     (void *) &init_flag);

        // Rebuild domain if necessary
        cout << "amrinit : Building new domain " << endl;
        fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);
        cout << "amrinit : Done building new domain " << endl << endl;

        if (new_domain != NULL)
        {
            // This is just for fun; remove when it gets annoying.
            // fclaw2d_domain_list_adapted(*domain, new_domain, SC_LP_STATISTICS);

            // Allocate memory for user data types (but they don't get set)
            allocate_user_data(new_domain);

            // Initialize new grids.  Assume that all ghost cells are filled in by qinit.
            fclaw2d_domain_iterate_adapted(*domain, new_domain,cb_domain_adapt,(void *) &init_flag);

            // Set some of the user data types.  Some of this is done in 'amr_set_base_level',
            // I should probably come up with a more general way to do this.
            set_domain_data(new_domain, gparms, amropts);
            set_domain_time(new_domain,t);

            // Physical BCs are needed in boundary level exchange
            // Assume only one block, since we are assuming mthbc

            int num = new_domain->num_blocks;
            for (int i = 0; i < num; i++)
            {
                fclaw2d_block_t *block = &new_domain->blocks[i];
                // This is kind of dumb for now, since block won't in general
                // have the same physical boundary conditions types.
                set_block_data(block,gparms->m_mthbc);
            }

            int new_level = level+1;
            // Upon initialization, we don't do any ghost cell exchanges, because we assume
            // that the initial conditions have set all internal ghost cells.
            bc_set_phys(new_domain,new_level,t);


            amrreset(*domain);
            fclaw2d_domain_destroy(*domain);

            *domain = new_domain;

            fclaw2d_domain_t *domain_partitioned =
                fclaw2d_domain_partition (*domain);
            if (domain_partitioned != NULL)
            {
                // TODO: allocate patch and block etc. memory for domain_partitioned

                // TODO: write a function to transfer values in parallel */

                /* then the old domain is no longer necessary */
                amrreset(*domain);
                fclaw2d_domain_destroy(*domain);
                *domain = domain_partitioned;

                /* internal clean up */
                fclaw2d_domain_complete(*domain);
            }
        }
        else
        {
            // exit loop;  we are done refining
            break;
        }
        cout << "amrinit : done with level " << endl << endl;;
    }
    cout << "Done with building initial grid structure " << endl;
}

void amrregrid(fclaw2d_domain_t **domain)
{

    fclaw2d_domain_data_t *ddata = get_domain_data(*domain);
    global_parms *gparms = ddata->parms;
    const amr_options_t *amropts = ddata->amropts;
    Real t = get_domain_time(*domain);

    bool init_flag = false;

    int minlevel = gparms->m_minlevel;
    int maxlevel = gparms->m_maxlevel;

    // Unlike the initial case, where we refine level by level, here, we only visit each tag
    // once and decide whether to refine or coarsen that patch.
    fclaw2d_domain_iterate_patches(*domain,
                                   (fclaw2d_patch_callback_t) cb_tag_patch,
                                   (void *) &init_flag);

    // Rebuild domain if necessary
    // Will return be NULL if no refining was done?
    cout << "amrregrid : Calling domain_adapt " << endl;
    fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);
    cout << "Done with domain adapt" << endl;

    if (new_domain != NULL)
    {
        // This is just for fun; remove when it gets annoying.
        // fclaw2d_domain_list_adapted(*domain, new_domain, SC_LP_STATISTICS);

        // Allocate memory for user data types (but they don't get set)
        allocate_user_data(new_domain);

        // Average or interpolate to new grids.
        fclaw2d_domain_iterate_adapted(*domain, new_domain,cb_domain_adapt,(void *) &init_flag);

        // Set some of the user data types.  Some of this is done in 'amr_set_base_level',
        // I should probably come up with a more general way to do this.
        set_domain_data(new_domain, gparms, amropts);
        set_domain_time(new_domain,t);

        // Physical BCs are needed in boundary level exchange
        // Assume only one block, since we are assuming mthbc
        int num = new_domain->num_blocks;
        for (int i = 0; i < num; i++)
        {
            fclaw2d_block_t *block = &new_domain->blocks[i];
            // This is kind of dumb for now, since block won't in general
            // have the same physical boundary conditions types.
            set_block_data(block,gparms->m_mthbc);
        }

        // Level stuff to make sure all
        for (int level = minlevel; level <= maxlevel; level++)
        {
            // Only do a level exchange;  Coarse and fine grid exchanges will happen
            // during time stepping as needed.
            bc_level_exchange(new_domain,level);
            bc_set_phys(new_domain,level,t);
        }

        amrreset(*domain);
        fclaw2d_domain_destroy(*domain);
        *domain = new_domain;

        fclaw2d_domain_t *domain_partitioned =
            fclaw2d_domain_partition (*domain);

        if (domain_partitioned != NULL)
        {
            // TODO: allocate patch and block etc. memory for domain_partitioned

            // TODO: write a function to transfer values in parallel */

            /* then the old domain is no longer necessary */
            amrreset(*domain);
            fclaw2d_domain_destroy(*domain);
            *domain = domain_partitioned;

            /* internal clean up */
            fclaw2d_domain_complete(*domain);
        }
    }
}


// -----------------------------------------------------------------
// Run - with or without subcycling
// -----------------------------------------------------------------
void amrrun(fclaw2d_domain_t *domain)
{

    int outstyle = 3;

    if (outstyle == 1)
    {
        explicit_step_fixed_output(domain);
    }
    else if (outstyle == 3)
    {
        int nstep = 20;  // Take 'nstep' steps
        int nplot = 1;   // Plot every 'nplot' steps
        explicit_step(domain,nstep,nplot);
    }
}

void explicit_step_fixed_output(fclaw2d_domain_t *domain)
{
    // Write out an initial time file
    int iframe = 0;
    amrout(domain,iframe);

    global_parms *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
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
            time_stepper.define(domain,gparms,ddata->amropts,t_curr);
            set_domain_time(domain,t_curr);

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
                printf("   WARNING : Took small time step which was %6.1f%% of desired dt.\n",
                       100.0*dt0/dt_level0);
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

            int regrid_step = 1;  // Will eventually read this in as a parameter.
            if (n_inner % regrid_step == 0)
            {
                // After some number of time steps, we probably need to regrid...
                amrregrid(&domain);
            }
        }

        // Output file at every outer loop iteration
        set_domain_time(domain,t_curr);
        iframe++;
        amrout(domain,iframe);
    }
}


void explicit_step(fclaw2d_domain_t *domain, int nstep_outer, int nstep_inner)
{
    // Write out an initial time file
    int iframe = 0;
    amrout(domain,iframe);

    global_parms *gparms = get_domain_parms(domain);
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    Real initial_dt = gparms->m_initial_dt;

    Real t0 = 0;

    // int nstep_outer = 2;  // Take this many steps in total
    // int nstep_inner = 1;   // Every 'nstep_inner', write out file.

    Real dt_level0 = initial_dt;
    Real t_curr = t0;
    set_domain_time(domain,t_curr);
    int n = 0;
    while (n < nstep_outer)
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

        // This also sets the time step on all finer levels.
        time_stepper.set_dt_minlevel(dt_minlevel);

        Real maxcfl_step = advance_all_levels(domain, &time_stepper);

        printf("Level %d step %5d : dt = %12.3e; maxcfl (step) = %8.3f; Final time = %12.4f\n",
               time_stepper.minlevel(),n+1,dt_minlevel,maxcfl_step, t_curr);

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
        n++;

        set_domain_time(domain,t_curr);

        // New time step, which should give a cfl close to the desired cfl.
        dt_level0 = dt_level0*gparms->m_desired_cfl/maxcfl_step;

        int regrid_step = 1;  // Will eventually read this in as a parameter.
        if (n % regrid_step == 0)
        {
            // After some number of time steps, we probably need to regrid...
            cout << "regridding at step " << n << endl;
            amrregrid(&domain);
        }

        if (n % nstep_inner == 0)
        {
            iframe++;
            amrout(domain,iframe);
        }
    }
}

void cb_amrout(fclaw2d_domain_t *domain,
               fclaw2d_patch_t *this_patch,
               int this_block_idx,
               int this_patch_idx,
               void *user)
{
    int iframe = *((int *) user);
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int num = this_block->num_patches_before + this_patch_idx + 1;
    int matlab_level = this_patch->level + 1;  // Matlab wants levels to start at 1.

    // Patch data is appended to fort.qXXXX
    ClawPatch *cp = get_clawpatch(this_patch);

    cp->write_patch_data(iframe, num, matlab_level);
}


void amrout(fclaw2d_domain_t *domain, int iframe)
{
    global_parms *gparms = get_domain_parms(domain);
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
