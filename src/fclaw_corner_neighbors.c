/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_corner_neighbors.h>

#include <fclaw_block.h>
#include <fclaw_ghost_fill.h>
#include <fclaw2d_map_query.h>
#include <fclaw_options.h>
#include <fclaw_global.h>
#include <fclaw_physical_bc.h>
#include <fclaw_patch.h>



/* This is used to determine neighbor patch relative level (finer, coarser or samesize)
   This enum is defined both here and in fclaw2d_face_neighbors.cpp.  Is that okay? */
enum
{
    COARSER_GRID = -1,
    SAMESIZE_GRID,
    FINER_GRID
};


static
int get_num_intersections(int dim,int intersects[], int faces[])
{
    int num_intersections = 0;
    int i;
    for (i = 0; i < dim; i++)
    {
        if (intersects[faces[i]])
        {
            num_intersections++;
        }
    }
    return num_intersections;
}

static
int find_face(int dim, int intersects[], int faces[])
{
    int i;
    for (i = 0; i < dim; i++)
    {
        if (intersects[faces[i]])
        {
            return faces[i];
        }
    }
    return -1;
}

/**
 * @brief Get the corner type for a patch
 * 
 * @param glob the global context
 * @param icorner the corner to get the cornertype on
 * @param intersects_bdry Size 4 for 2D, 6 for 3D. True if patch intersects a boundary on a face
 * @param intersects_block Size 4 for 2D, 6 for 3D. True if patch intersects a block boundary on a face
 * @param interior_corner Returns true if the corner is interior to a block
 * @param ftransform The ftransform object to fill in
 *                      This will fill in the is_block_corner, is_block_edge, is_block_face, block_iface, block_iedge, and block_iface fields.
 * @param ftransform_fine The ftransform object to fill in
 *                      This will fill in the is_block_corner, is_block_edge, and is_block_face fields.
 *                      block_iface and block_iedge will be initialized to -1.
 */
static
void get_corner_type(fclaw_global_t* glob,
                     int icorner,
                     int intersects_bdry[],
                     int intersects_block[],
                     int *is_interior_in_domain,
                     fclaw_patch_transform_data_t* tdata,
                     fclaw_patch_transform_data_t* tdata_finegrid)
{
    fclaw_domain_t *domain = glob->domain;

    tdata->icorner = icorner;
    tdata->iedge = -1;
    tdata->iface = -1;

    // p4est has tons of lookup table like this, can be exposed similarly
    int corner_faces[domain->refine_dim];
    fclaw_domain_corner_faces(domain, icorner, corner_faces);

    /* Both faces are at a physical boundary */
    int num_phys_faces = get_num_intersections(domain->refine_dim, 
                                               intersects_bdry,
                                               corner_faces);

    /* Either a corner is at a block boundary (but not a physical boundary),
       or internal to a block.  L-shaped domains are excluded for now
       (i.e. no reentrant corners). */
    *is_interior_in_domain = num_phys_faces == 0;

    int num_block_faces = get_num_intersections(domain->refine_dim, 
                                                intersects_block,
                                                corner_faces);
    /* Both faces are at a block boundary, physical or not */
    tdata->is_block_corner = num_block_faces == domain->refine_dim;

    tdata->is_block_face = num_block_faces == 1;
    tdata->block_iface = -1;
    if (num_block_faces == 1)
    {
        tdata->block_iface = find_face(domain->refine_dim, intersects_block, corner_faces);
    }

    tdata_finegrid->is_block_corner = tdata->is_block_corner;
    tdata_finegrid->is_block_edge = tdata->is_block_edge;
    tdata_finegrid->is_block_face = tdata->is_block_face;

    tdata_finegrid->block_iface = -1;
    tdata_finegrid->block_iedge = -1;
}

/* --------------------------------------------------------
   Four cases to consider.   The 'has_corner_neighbor'
   value is returned from p4est.  The assumption going
   into this routine is that we have found a valid
   interior corner (not a corner on a physical boundary).
   The corner than satisfies one of the following four
   cases.

   Case No. | has_corner_neighbor  |  block_boundary_type
   --------------------------------------------------------
      1     |       T              |    interior
      2     |       T              |    block_corner
      3     |       T              |    block_edge
      4     |       T              |    block_face
      5     |       F              |    interior
      6     |       F              |    block_corner
      7     |       F              |    block_edge
      8     |       F              |    block_face

    Case 1 : Interior corner 
    Case 2 : Corner neighbor across block corner
    Case 3 : Corner neighbor across block edge
    Case 4 : Corner neighbor across block face
    Case 5-8 : Invalid, not curerntly supported
   ------------------------------------------------------ */

static
void get_corner_neighbor_3d(fclaw_global_t *glob,
                            int this_block_idx,
                            int this_patch_idx,
                            fclaw_patch_t* this_patch,
                            int icorner,
                            int block_iface,
                            int *corner_block_idx,
                            fclaw_patch_t** corner_patch,
                            int *rcornerno,
                            int **ref_flag_ptr,
                            int *block_corner_count,
                            fclaw_patch_transform_data_t* ftransform,
                            fclaw_patch_transform_data_t* ftransform_finegrid)
{
    fclaw_domain_t *domain = glob->domain;
    /* See what p4est thinks we have for corners, and consider four cases */
    int rproc_corner;
    int corner_patch_idx;
    fclaw_patch_relation_t neighbor_type;

    /* Note : Pillowsphere case does not return a block corner neighbor */
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);
    int has_corner_neighbor =
        fclaw_patch_corner_neighbors(domain,
                                       this_block_idx,
                                       this_patch_idx,
                                       icorner,
                                       &rproc_corner,
                                       corner_block_idx,
                                       &corner_patch_idx,
                                       rcornerno,
                                       &neighbor_type);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);    

    *block_corner_count = 0;  /* Assume we are not at a block corner */
    if (has_corner_neighbor)
    {
        /* Cases 1-4: */
        *block_corner_count = 4;  /* assume four for now */

        /* TODO: use interior for now, this needs to be changed when we want to use something other than brick */
        fclaw_patch_transform_blockface_intra (glob, ftransform->transform);
        fclaw_patch_transform_blockface_intra (glob, ftransform_finegrid->transform);
    }
    else
    {
        /* Cases 5-8 : not supported */
        /* We do not return valid transformation objects! */
        *ref_flag_ptr = NULL;
        *corner_patch = NULL;
        return;
    }

    /* ---------------------------------------------------------------------
       We have a valid neighbor and possibly a transform. We just now need
       to get a pointer to the neighbor patch (which may be a parallel patch)
       and the relative level (-1,0,1).
       --------------------------------------------------------------------- */

    if (domain->mpirank != rproc_corner)
    {
        *corner_patch = &domain->ghost_patches[corner_patch_idx];
    }
    else
    {
        fclaw_block_t *neighbor_block = &domain->blocks[*corner_block_idx];
        *corner_patch = &neighbor_block->patches[corner_patch_idx];
    }

    if (neighbor_type == FCLAW_PATCH_HALFSIZE)
    {
        **ref_flag_ptr = 1;
    }
    else if (neighbor_type == FCLAW_PATCH_SAMESIZE)
    {
        **ref_flag_ptr = 0;
    }
    else /* FCLAW2D_PATCH_DOUBLESIZE */
    {
        **ref_flag_ptr = -1;
    }
}
/* --------------------------------------------------------
   Four cases to consider.   The 'has_corner_neighbor'
   value is returned from p4est.  The assumption going
   into this routine is that we have found a valid
   interior corner (not a corner on a physical boundary).
   The corner than satisfies one of the following four
   cases.

   Case No. | has_corner_neighbor  |  is_block_corner
   --------------------------------------------------------
      1     |       T              |        T
      2     |       F              |        F
      3     |       T              |        F
      4     |       F              |        T

    Case 1 : In this case, 4 or more patches meet at a block
             corner. No transforms are yet available, so we 
             assume that at block corners, the patches all 
             have the same orientation.
    Case 2 : Corner is at a hanging node and has no valid
             adjacent corner.
    Case 3 : Corner is either interior to a block, or on a
             block edge.  In each case, the transform is
             well-defined.
    Case 4 : Either 3 patches meet at a corner, in which
             case we don't have an adjacent corner, or we are
             on a pillow grid, in which case we have a 
             corner, but one which we nonetheless treat
             as a special case.
   ------------------------------------------------------ */

static
void get_corner_neighbor_2d(fclaw_global_t *glob,
                            int this_block_idx,
                            int this_patch_idx,
                            fclaw_patch_t* this_patch,
                            int icorner,
                            int block_iface,
                            int *corner_block_idx,
                            fclaw_patch_t** corner_patch,
                            int *rcornerno,
                            int **ref_flag_ptr,
                            int *block_corner_count,
                            fclaw_patch_transform_data_t* ftransform,
                            fclaw_patch_transform_data_t* ftransform_finegrid)
{
    fclaw_domain_t *domain = glob->domain;
    /* See what p4est thinks we have for corners, and consider four cases */
    int rproc_corner;
    int corner_patch_idx;
    fclaw_patch_relation_t neighbor_type;

    /* Note : Pillowsphere case does not return a block corner neighbor */
    int ispillowsphere = 0;
    if(domain->refine_dim == 2)
    {
        //TODO 3d
        fclaw2d_map_pillowsphere(glob);
    }

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);
    int has_corner_neighbor =
        fclaw_patch_corner_neighbors(domain,
                                       this_block_idx,
                                       this_patch_idx,
                                       icorner,
                                       &rproc_corner,
                                       corner_block_idx,
                                       &corner_patch_idx,
                                       rcornerno,
                                       &neighbor_type);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);    

    *block_corner_count = 0;  /* Assume we are not at a block corner */
    if (has_corner_neighbor && ftransform->is_block_corner)
    {
        /* Case 1 : 4 or more patches meet at a block corner.
        This case does NOT include the pillowgrid.   */
        *block_corner_count = 4;  /* assume four for now */

        /* No block corner transforms yet, so we use the 
        interior 'default' transforms. */
        fclaw_patch_transform_blockface_intra (glob, ftransform->transform);
        fclaw_patch_transform_blockface_intra
            (glob, ftransform_finegrid->transform);
    }
    else if (!has_corner_neighbor && !ftransform->is_block_corner)
    {
        /* Case 2 : 'icorner' is a hanging node */
        /* We do not return valid transformation objects! */
        *ref_flag_ptr = NULL;
        *corner_patch = NULL;
        return;
    }
    else if (has_corner_neighbor && !ftransform->is_block_corner)
    {
        /* Case 3 : 'icorner' is an interior corner, at a block edge,
         or we are on a periodic block.  Need to return a valid
         transform in 'ftransform' */
        /* block_iface is the block number at the of the neighbor? */
        if (block_iface >= 0)
        {
            /* The corner is on a block edge (but is not a block corner).
               Compute a transform between blocks. First, get the
               remote face number.  The remote face number encodes the
               orientation, so we have 0 <= rfaceno < 8 */
            int rfaceno;
            int rproc[4]; // overallocate for 3d
            int rpatchno[4];
            int rblockno;  /* Should equal *corner_block_idx, above. */
            fclaw_patch_face_neighbors(domain,
                                         this_block_idx,
                                         this_patch_idx,
                                         block_iface,
                                         rproc,
                                         &rblockno,
                                         rpatchno,
                                         &rfaceno);

            FCLAW_ASSERT(rblockno == *corner_block_idx);

            /* Get encoding of transforming a neighbor coordinate across a face */
            fclaw_patch_transform_blockface (glob, block_iface, rfaceno, ftransform->transform);

            /* Get transform needed to swap parallel ghost patch with fine
               grid on-proc patch.  This is done so that averaging and
               interpolation routines can be re-used. */
            int iface1 = block_iface;
            int rface1 = rfaceno;
            fclaw_patch_face_swap(domain->refine_dim, &iface1, &rface1);
            fclaw_patch_transform_blockface(glob, iface1, rface1,
                                              ftransform_finegrid->transform);

            ftransform_finegrid->block_iface = iface1;
        }
        else if (this_block_idx == *corner_block_idx)
        {
            /* Both patches are in the same block, so we set the transform to
               a default transform.  This could be the case for periodic boundaries. */
            *block_corner_count = 4;  /* assume four for now */
            fclaw_patch_transform_blockface_intra (glob, ftransform->transform);
            fclaw_patch_transform_blockface_intra
                (glob, ftransform_finegrid->transform);

        }
        else
        {
            fclaw_errorf("WARNING : this_block_idx = %d\n",this_block_idx);
            fclaw_errorf("WARNING : corner_block_idx = %d\n",*corner_block_idx);
            fclaw_abortf("get_corner_neighbors " \
                            "(fclaw2d_corner_neighbors.c : " \
                            "We should not be here\n");
        }
    }
    else if (!has_corner_neighbor && ftransform->is_block_corner)
    {
        /* Case 4 : Pillow sphere case or cubed sphere  */
        if (!ispillowsphere)
        {
            *block_corner_count = 3;
            /* Exactly 3 patches meet at a corner, e.g. the cubed sphere.
               In this case, 'this_patch' has no corner-adjacent only
               neighbors, and so there is nothing to do. */
            *ref_flag_ptr = NULL;
            *corner_patch = NULL;
            return;
        }
        else
        {
            *block_corner_count = 2;
            has_corner_neighbor = 1;
            int rpatchno[4]; // overallocate for 3d
            int rproc[4];
            int rfaceno;

            /* Use only faces 0 or 1 to get block data. */
            int iface = icorner % 2;
            neighbor_type =
                fclaw_patch_face_neighbors(domain,
                                             this_block_idx,
                                             this_patch_idx,
                                             iface,
                                             rproc,
                                             corner_block_idx,
                                             rpatchno,
                                             &rfaceno);

            int igrid;
            if (neighbor_type == FCLAW_PATCH_HALFSIZE)
            {
                /* igrid = 0 at corners 0,1 and (R-1) at corners 2,3,
                   where R = refinement factor */
                igrid = (icorner/2)*(2 - 1);
            }
            else
            {
                /* Same size or double size grids. */
                igrid = 0;
            }

            *rcornerno = icorner;    /* This wasn't being set! */
            corner_patch_idx = rpatchno[igrid];
            rproc_corner = rproc[igrid];
        }
    }

    /* ---------------------------------------------------------------------
       We have a valid neighbor and possibly a transform. We just now need
       to get a pointer to the neighbor patch (which may be a parallel patch)
       and the relative level (-1,0,1).
       --------------------------------------------------------------------- */

    if (domain->mpirank != rproc_corner)
    {
        *corner_patch = &domain->ghost_patches[corner_patch_idx];
    }
    else
    {
        fclaw_block_t *neighbor_block = &domain->blocks[*corner_block_idx];
        *corner_patch = &neighbor_block->patches[corner_patch_idx];
    }

    if (neighbor_type == FCLAW_PATCH_HALFSIZE)
    {
        **ref_flag_ptr = 1;
    }
    else if (neighbor_type == FCLAW_PATCH_SAMESIZE)
    {
        **ref_flag_ptr = 0;
    }
    else /* FCLAW2D_PATCH_DOUBLESIZE */
    {
        **ref_flag_ptr = -1;
    }
}



void cb_corner_fill(fclaw_domain_t *domain,
                    fclaw_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    void *user)
{
    const int num_faces = fclaw_domain_num_faces(domain);
    const int num_corners = fclaw_domain_num_corners(domain);

    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user; 

    fclaw_exchange_info_t *filltype = (fclaw_exchange_info_t*) s->user;
    int time_interp = filltype->time_interp;
    int is_coarse = filltype->grid_type == FCLAW_IS_COARSE;
    int is_fine = filltype->grid_type == FCLAW_IS_FINE;

    int read_parallel_patches = filltype->read_parallel_patches;

    int copy_from_neighbor = filltype->exchange_type == FCLAW_COPY;
    int average_from_neighbor = filltype->exchange_type == FCLAW_AVERAGE;
    int interpolate_to_neighbor = filltype->exchange_type == FCLAW_INTERPOLATE;


    int intersects_bdry[num_faces];
    fclaw_physical_get_bc(s->glob,this_block_idx,this_patch_idx,
                            intersects_bdry);

    int intersects_block[num_faces];
    fclaw_block_get_block_boundary(s->glob, this_patch, intersects_block);

    /* Transform data needed at multi-block boundaries */
    fclaw_patch_transform_data_t transform_data;

    transform_data.glob = s->glob;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.

    fclaw_patch_transform_init_data(s->glob,this_patch,
                                      this_block_idx,
                                      this_patch_idx,
                                      &transform_data);


    fclaw_patch_transform_data_t transform_data_finegrid;

    transform_data_finegrid.glob = s->glob;
    transform_data_finegrid.based = 1;   // cell-centered data in this routine.

    fclaw_patch_transform_init_data(s->glob,this_patch,
                                      this_block_idx,
                                      this_patch_idx,
                                      &transform_data_finegrid);


    for (int icorner = 0; icorner < num_corners; icorner++)
    {
        int block_corner_count = 0;
        // get corner type and initialize is_block_* values in transform_data and transform_data_finegrid
        // block_i* values will be initialized for transform_data
        // block_i* values will be initialized to -1 for transform_data_finegrid
        int is_interior_in_domain;
        get_corner_type(s->glob,icorner,
                        intersects_bdry,
                        intersects_block,
                        &is_interior_in_domain,
                        &transform_data,
                        &transform_data_finegrid);

        /* Sets block_corner_count to 0 */
        fclaw_patch_set_block_corner_count(s->glob, this_patch,
                                             icorner,block_corner_count);

        if (is_interior_in_domain)
        {
            /* Is an interior patch corner;  may also be a block corner */

            int corner_block_idx;
            int neighbor_level;
            int *ref_flag_ptr = &neighbor_level;
            fclaw_patch_t *corner_patch;
            int rcornerno;

            int  block_iface = transform_data.block_iface;

            corner_block_idx = -1;
            get_corner_neighbor_2d(s->glob,
                                this_block_idx,
                                this_patch_idx,
                                this_patch,
                                icorner,
                                block_iface,
                                &corner_block_idx,
                                &corner_patch,
                                &rcornerno,
                                &ref_flag_ptr,
                                &block_corner_count,
                                &transform_data,
                                &transform_data_finegrid);

            /* This sets value in block_corner_count_array */
            fclaw_patch_set_block_corner_count(s->glob, this_patch,
                                                 icorner,block_corner_count);


            /* Needed for switching the context */
            transform_data_finegrid.icorner = rcornerno;
            transform_data_finegrid.this_patch = corner_patch;
            transform_data_finegrid.neighbor_patch = this_patch;
            
            if (ref_flag_ptr == NULL)
            {
                /* No corner neighbor.  Either :
                   -- Hanging node
                   -- Cubed sphere
                */
                continue;
            }

            int remote_neighbor = fclaw_patch_is_ghost(corner_patch);
            if (is_coarse && ((read_parallel_patches && remote_neighbor) || !remote_neighbor))
            {
                transform_data.neighbor_patch = corner_patch;
                if (neighbor_level == FINER_GRID)
                {
                    fclaw_patch_t* coarse_patch = this_patch;
                    fclaw_patch_t* fine_patch = corner_patch;
                    int coarse_blockno = this_block_idx;
                    int fine_blockno = corner_block_idx;
                    if (interpolate_to_neighbor && !remote_neighbor)
                    {
                        /* No need to interpolate to remote ghost patches. */
                        fclaw_patch_interpolate_corner(s->glob,
                                                       coarse_patch,
                                                       fine_patch,
                                                       coarse_blockno,
                                                       fine_blockno,
                                                       transform_data.is_block_corner,
                                                       icorner,time_interp,
                                                       &transform_data);
                    }
                    else if (average_from_neighbor)
                    {
                        /* Average even if neighbor is a remote neighbor */
                        fclaw_patch_t* coarse_patch = this_patch;
                        fclaw_patch_t* fine_patch = corner_patch;
                        fclaw_patch_average_corner(s->glob,
                                                   coarse_patch,
                                                   fine_patch,
                                                   coarse_blockno,
                                                   fine_blockno,
                                                   transform_data.is_block_corner,
                                                   icorner,time_interp,
                                                   &transform_data);                        
                    }
                }
                else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                {
                    fclaw_patch_copy_corner(s->glob,
                                            this_patch,
                                            corner_patch,
                                            this_block_idx,
                                            corner_block_idx,
                                            transform_data.is_block_corner,
                                            icorner, time_interp,
                                            &transform_data);
                }

            }  /* End of non-parallel patch case */
            else if (is_fine && neighbor_level == COARSER_GRID &&
                     remote_neighbor && read_parallel_patches)
            {
                /* The coarse grid is now the remote patch;  swap contexts and 
                call same routines above, but with remote patch as the "coarse" 
                grid */
                
                fclaw_patch_t* coarse_patch = corner_patch;
                fclaw_patch_t* fine_patch = this_patch;
                int coarse_blockno = corner_block_idx;
                int fine_blockno = this_patch_idx;

                if (interpolate_to_neighbor)
                {
                    /* Interpolate from remote coarse grid patch (coarse grid) to
                       local fine grid patch.  We do not need to average to the 
                       remote patch corners unless corners are used in the 
                       interpolation stencil. */
                    int coarse_icorner = transform_data_finegrid.icorner;

                    fclaw_patch_interpolate_corner(s->glob,
                                                   coarse_patch,
                                                   fine_patch,
                                                   coarse_blockno,
                                                   fine_blockno,
                                                   transform_data.is_block_corner,
                                                   coarse_icorner,time_interp,
                                                   &transform_data_finegrid);

                }
            } /* End of parallel case */
        }  /* End of 'interior_corner' */
    }  /* End of icorner loop */
}
