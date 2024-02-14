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
#include <fclaw_map_query.h>
#include <fclaw_options.h>
#include <fclaw_global.h>
#include <fclaw_physical_bc.h>
#include <fclaw_patch.h>



static
int get_num_intersections(int dim,int intersects[], int faces[])
{
    int num_intersections = 0;
    for (int i = 0; i < dim; i++)
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
    for (int i = 0; i < dim; i++)
    {
        if (intersects[faces[i]])
        {
            return faces[i];
        }
    }
    return -1;
}

static
int find_edge(int corner, int intersects[], int faces[])
{
    int non_intersecting_face = -1;
    for (int i = 0; i < 3; i++)
    {
        if (!intersects[faces[i]])
        {
            non_intersecting_face = faces[i];
            break;
        }
    }
    FCLAW_ASSERT(non_intersecting_face >= 0);
    int edge_axis = non_intersecting_face / 2;
    int face_0;
    int face_1;
    switch(edge_axis)
    {
        case 0:
            face_0 = faces[1];
            face_1 = faces[2];
            break;
        case 1:
            face_0 = faces[0];
            face_1 = faces[2];
            break;
        case 2:
            face_0 = faces[0];
            face_1 = faces[1];
            break;
    }
    FCLAW_ASSERT(face_0 != non_intersecting_face);
    FCLAW_ASSERT(face_1 != non_intersecting_face);

    int upper_0 = face_0 % 2;
    int upper_1 = face_1 % 2;

    int edge = 4*edge_axis + 2*upper_1 + upper_0;

    FCLAW_ASSERT(edge >= 0 && edge < 12);

    return edge;
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

    if(domain->refine_dim == 2)
    {
        tdata->is_block_edge = 0;
        tdata->block_iedge = -1;
    }
    else
    {
        tdata->is_block_edge = num_block_faces == 2;
        tdata->block_iedge = -1;
        if(tdata->is_block_edge)
        {
        tdata->block_iedge = find_edge(icorner, intersects_block, corner_faces);
        }
    }

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
      1     |       T              |    is_block_corner
      2     |       T              |    is_block_edge
      3     |       T              |    is_block_face
      4     |       T              |    interior
      2     |       F              |        F
      3     |       T              |        F
      4     |       F              |        T

    Case 1-2 : In this case, 4 or more patches meet at a
             corner. No transforms are yet available, so we 
             assume that at block corners, the patches all 
             have the same orientation.
    Case 3 : Corner is either is on block face.
             The transform is well-defined.
    Case 4 : Corner is either interior to a block.
             The transform is well-defined.
    Case 2 : Corner is at a hanging node and has no valid
             adjacent corner.
    Case 4 : Either 3 patches meet at a corner, in which
             case we don't have an adjacent corner, or we are
             on a pillow grid, in which case we have a 
             corner, but one which we nonetheless treat
             as a special case.
   ------------------------------------------------------ */

static
void get_corner_neighbor(int indirect, 
                         fclaw_global_t *glob,
                         int icorner,
                         int *is_valid_neighbor,
                         int *block_corner_count,
                         fclaw_patch_transform_data_t* tdata,
                         fclaw_patch_transform_data_t* tdata_fine)
{
    fclaw_domain_t *domain = glob->domain;
    /* assume neighbor is valid for now*/
    *is_valid_neighbor = 1;
    /* See what p4est thinks we have for corners, and consider four cases */
    int rproc_corner;
    
    /* Note : Pillowsphere case does not return a block corner neighbor */
    int ispillowsphere = 0;
    if(domain->refine_dim == 2)
    {
        ispillowsphere =  FCLAW_MAP_IS_PILLOWSPHERE(&glob->cont);
    }

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);
    int has_corner_neighbor = 0;
    if(indirect)
    {
        tdata->neighbor_type =
            fclaw_domain_indirect_corner_neighbor(domain,
                                                  domain->indirect,
                                                  tdata->this_patchno,
                                                  icorner,
                                                  &rproc_corner,
                                                  &tdata->neighbor_blockno,
                                                  &tdata->neighbor_patchno,
                                                  &tdata_fine->icorner);
        has_corner_neighbor = tdata->neighbor_type != FCLAW_PATCH_BOUNDARY;
    }
    else
    {
        has_corner_neighbor =
            fclaw_patch_corner_neighbors(domain,
                                         tdata->this_blockno,
                                         tdata->this_patchno,
                                         icorner,
                                         &rproc_corner,
                                         &tdata->neighbor_blockno,
                                         &tdata->neighbor_patchno,
                                         &tdata_fine->icorner,
                                         &tdata->neighbor_type);
    }

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);    

    *block_corner_count = 0;  /* Assume we are not at a block corner */
    if (has_corner_neighbor && tdata->is_block_corner)
    {
        /* Case 1 : 4 or more patches meet at a block corner.
        This case does NOT include the pillowgrid.   */
        *block_corner_count = 4;  /* assume four for now */

        /* No block corner transforms yet, so we use the 
        interior 'default' transforms. */
        fclaw_patch_transform_blockface_intra (glob, tdata->transform);
        fclaw_patch_transform_blockface_intra
            (glob, tdata_fine->transform);
    }
    else if (has_corner_neighbor && tdata->is_block_edge)
    {
        FCLAW_ASSERT(domain->refine_dim == 3);
        /* The corner is on a block face (but is not a block corner).
           Compute a transform between blocks. First, get the
           remote face number.  The remote face number encodes the
           orientation, so we have 0 <= rfaceno < 8 */

        //int rfaceno;
        //int rproc[2];
        //int rpatchno[2];
        //int rblockno;  /* Should equal *corner_block_idx, above. */
        //fclaw_patch_edge_neighbors(domain,
        //                             tdata->this_blockno,
        //                             tdata->this_patchno,
        //                             tdata->block_iface,
        //                             rproc,
        //                             &rblockno,
        //                             rpatchno,
        //                             &rfaceno);

        //FCLAW_ASSERT(rblockno == tdata->neighbor_blockno);

        ///* Get encoding of transforming a neighbor coordinate across a face */
        //fclaw_patch_transform_blockface (glob, tdata->block_iface, rfaceno, tdata->transform);
        fclaw_patch_transform_blockface_intra(glob, tdata->transform);

        /* Get transform needed to swap parallel ghost patch with fine
           grid on-proc patch.  This is done so that averaging and
           interpolation routines can be re-used. */
        //int iface1 = tdata->block_iface;
        //int rface1 = rfaceno;
        //fclaw_patch_face_swap(domain->refine_dim, &iface1, &rface1);
        //TODO edge transforms
        fclaw_patch_transform_blockface_intra(glob, tdata_fine->transform);

        //TODO update this for rotation
        tdata_fine->block_iedge = tdata->block_iedge ^ 3;
    }
    else if (has_corner_neighbor && tdata->is_block_face)
    {
        /* The corner is on a block face (but is not a block corner).
           Compute a transform between blocks. First, get the
           remote face number.  The remote face number encodes the
           orientation, so we have 0 <= rfaceno < 8 */

        int rfaceno;
        int rproc[4]; // overallocate for 3d
        int rpatchno[4];
        int rblockno;  /* Should equal *corner_block_idx, above. */
        if(indirect)
        {
            fclaw_domain_indirect_face_neighbors(domain, 
                                                domain->indirect, 
                                                tdata->this_patchno, 
                                                tdata->block_iface, 
                                                rproc, 
                                                &rblockno, 
                                                rpatchno,
                                                &rfaceno);
        }
        else 
        {
            fclaw_patch_face_neighbors(domain,
                                       tdata->this_blockno,
                                       tdata->this_patchno,
                                       tdata->block_iface,
                                       rproc,
                                       &rblockno,
                                       rpatchno,
                                       &rfaceno);
            FCLAW_ASSERT(rblockno == tdata->neighbor_blockno);
        }


        /* Get encoding of transforming a neighbor coordinate across a face */
        fclaw_patch_transform_blockface (glob, tdata->block_iface, rfaceno, tdata->transform);

        /* Get transform needed to swap parallel ghost patch with fine
           grid on-proc patch.  This is done so that averaging and
           interpolation routines can be re-used. */
        int iface1 = tdata->block_iface;
        int rface1 = rfaceno;
        fclaw_patch_face_swap(domain->refine_dim, &iface1, &rface1);
        fclaw_patch_transform_blockface(glob, iface1, rface1,
                                        tdata_fine->transform);
        
        tdata_fine->block_iface = iface1;
    }
    else if (has_corner_neighbor && !tdata->is_block_corner)
    {
        /* Case 4 : 'icorner' is an interior corner,
         or we are on a periodic single block.  Need to return a valid
         transform in 'ftransform' */
        if (tdata->this_blockno == tdata->neighbor_blockno)
        {
            /* Both patches are in the same block, so we set the transform to
               a default transform.  This could be the case for periodic boundaries. */
            *block_corner_count = 4;  /* assume four for now */
            fclaw_patch_transform_blockface_intra (glob, tdata->transform);
            fclaw_patch_transform_blockface_intra (glob, tdata_fine->transform);

        }
        else
        {
            fclaw_errorf("WARNING : this_block_idx = %d\n",tdata->this_blockno);
            fclaw_errorf("WARNING : corner_block_idx = %d\n",tdata->neighbor_blockno);
            fclaw_abortf("get_corner_neighbors " \
                                    "(fclaw2d_corner_neighbors.c : " \
                                    "We should not be here\n");
        }
    }
    else if (!has_corner_neighbor && tdata->is_block_corner)
    {   
        /* not needed for indirect */
        if(indirect)
        {
            *block_corner_count = 0;
            *is_valid_neighbor = 0;
            return;
        }

        /* Case 4 : In 2d: Pillow sphere case or cubed sphere
           In 3D: not yet supported */
        if(domain->refine_dim == 3)
        {
            *block_corner_count = 0;
            *is_valid_neighbor = 0;
            return;
        }

        if (!ispillowsphere)
        {
            *block_corner_count = 3;
            /* Exactly 3 patches meet at a corner, e.g. the cubed sphere.
               In this case, 'this_patch' has no corner-adjacent only
               neighbors, and so there is nothing to do. */
            *is_valid_neighbor = 0;
            return;
        }
        else
        {
            /* return face neighbor as corner neighbor */
            *block_corner_count = 2;
            has_corner_neighbor = 1;
            int rpatchno[4]; // overallocate for 3d
            int rproc[4];
            int rfaceno;

            /* Use only faces 0 or 1 to get block data. */
            int iface = icorner % 2;
            tdata->neighbor_type =
                fclaw_patch_face_neighbors(domain,
                                             tdata->this_blockno,
                                             tdata->this_patchno,
                                             iface,
                                             rproc,
                                             &tdata->neighbor_blockno,
                                             rpatchno,
                                             &rfaceno);

            int igrid;
            if (tdata->neighbor_type == FCLAW_PATCH_HALFSIZE)
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

            tdata_fine->icorner = icorner;    /* This wasn't being set! */
            tdata->neighbor_patchno = rpatchno[igrid];
            rproc_corner = rproc[igrid];
        }
    }
    else if (!has_corner_neighbor)
    {
        /* 'icorner' is a hanging node or not yet supported */
        /* We do not return valid transformation objects! */
        *is_valid_neighbor = 0;
        return;
    }

    /* ---------------------------------------------------------------------
       We have a valid neighbor and possibly a transform. We just now need
       to get a pointer to the neighbor patch (which may be a parallel patch)
       and the relative level (-1,0,1).
       --------------------------------------------------------------------- */

    if (domain->mpirank != rproc_corner)
    {
        tdata->neighbor_patch = &domain->ghost_patches[tdata->neighbor_patchno];
    }
    else
    {
        fclaw_block_t *neighbor_block = &domain->blocks[tdata->neighbor_blockno];
        tdata->neighbor_patch = &neighbor_block->patches[tdata->neighbor_patchno];
    }

    

    // set neighbor information in tdata_finegrid
    tdata_fine->this_patch = tdata->neighbor_patch;
    tdata_fine->this_blockno = tdata->neighbor_blockno;
    tdata_fine->this_patchno = tdata->neighbor_patchno;

    if (tdata->neighbor_type == FCLAW_PATCH_HALFSIZE)
    {
        tdata_fine->neighbor_type = FCLAW_PATCH_DOUBLESIZE;
    }
    else if (tdata->neighbor_type == FCLAW_PATCH_SAMESIZE)
    {
        tdata_fine->neighbor_type = FCLAW_PATCH_SAMESIZE;
    }
    else /* FCLAW2D_PATCH_DOUBLESIZE */
    {
        tdata_fine->neighbor_type = FCLAW_PATCH_HALFSIZE;
    }
}



void fclaw_corner_fill_cb(fclaw_domain_t *domain,
                          fclaw_patch_t *this_patch,
                          int this_blockno,
                          int this_patchno,
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
    fclaw_physical_get_bc(s->glob,this_blockno,this_patchno,
                            intersects_bdry);

    int intersects_block[num_faces];
    fclaw_block_get_block_boundary(s->glob, this_patch, intersects_block);

    /* Transform data needed at multi-block boundaries */
    fclaw_patch_transform_data_t tdata;

    tdata.glob = s->glob;
    tdata.based = 1;   // cell-centered data in this routine.
    tdata.this_patch = this_patch;
    tdata.this_blockno = this_blockno;
    tdata.this_patchno = this_patchno;
    tdata.neighbor_patch = NULL;  // gets filled in below.

    fclaw_patch_transform_init_data(s->glob,this_patch,
                                      this_blockno,
                                      this_patchno,
                                      &tdata);


    fclaw_patch_transform_data_t tdata_fine;

    tdata_fine.glob = s->glob;
    tdata_fine.neighbor_patch = this_patch;
    tdata_fine.neighbor_blockno = this_blockno;
    tdata_fine.neighbor_patchno = this_patchno;
    tdata_fine.based = 1;   // cell-centered data in this routine.

    fclaw_patch_transform_init_data(s->glob,this_patch,
                                      this_blockno,
                                      this_patchno,
                                      &tdata_fine);


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
                        &tdata,
                        &tdata_fine);

        /* Sets block_corner_count to 0 */
        fclaw_patch_set_block_corner_count(s->glob, this_patch,
                                             icorner,block_corner_count);

        if (is_interior_in_domain)
        {
            /* Is an interior patch corner;  may also be a block corner */

            int is_valid_neighbor;

            int indirect = 0;
            get_corner_neighbor(indirect,
                                s->glob,
                                icorner,
                                &is_valid_neighbor,
                                &block_corner_count,
                                &tdata,
                                &tdata_fine);

            /* This sets value in block_corner_count_array */
            fclaw_patch_set_block_corner_count(s->glob, this_patch,
                                                 icorner,block_corner_count);
            

            if (!is_valid_neighbor)
            {
                /* No corner neighbor.  Either :
                   -- Hanging node
                   -- Cubed sphere
                */
                continue;
            }

            int remote_neighbor = fclaw_patch_is_ghost(tdata.neighbor_patch);
            if (is_coarse)
            {
                if (tdata.neighbor_type == FCLAW_PATCH_HALFSIZE)
                {
                    if (interpolate_to_neighbor && !remote_neighbor)
                    {
                        /* No need to interpolate to remote ghost patches. */
                        fclaw_patch_interpolate_corner(s->glob,
                                                       tdata.this_patch,
                                                       tdata.neighbor_patch,
                                                       tdata.this_blockno,
                                                       tdata.neighbor_blockno,
                                                       tdata.is_block_corner,
                                                       tdata.icorner,
                                                       time_interp,
                                                       &tdata);
                    }
                    else if (average_from_neighbor)
                    {
                        /* Average even if neighbor is a remote neighbor */
                        fclaw_patch_average_corner(s->glob,
                                                   tdata.this_patch,
                                                   tdata.neighbor_patch,
                                                   tdata.this_blockno,
                                                   tdata.neighbor_blockno,
                                                   tdata.is_block_corner,
                                                   tdata.icorner,
                                                   time_interp,
                                                   &tdata);                        
                    }
                }
                else if (tdata.neighbor_type == FCLAW_PATCH_SAMESIZE && copy_from_neighbor)
                {
                    fclaw_patch_copy_corner(s->glob,
                                            tdata.this_patch,
                                            tdata.neighbor_patch,
                                            tdata.this_blockno,
                                            tdata.neighbor_blockno,
                                            tdata.is_block_corner,
                                            tdata.icorner, 
                                            time_interp,
                                            &tdata);
                    
                    /* We also need to copy _to_ the remote neighbor; switch contexts, but
                       use ClawPatches that are only in scope here, to avoid
                       conflicts with above uses of the same variables. This is needed
                       in case we want to interpolate to adjacent corners on fine grids.*/
                    if (remote_neighbor)
                    {
                        fclaw_patch_copy_corner(s->glob,
                                                tdata_fine.this_patch,
                                                tdata_fine.neighbor_patch,
                                                tdata_fine.this_blockno,
                                                tdata_fine.neighbor_blockno,
                                                tdata_fine.is_block_corner,
                                                tdata_fine.icorner, 
                                                time_interp,
                                                &tdata_fine);                            
                    }
                }

            }  /* End of non-parallel patch case */
            else if (is_fine && tdata.neighbor_type == FCLAW_PATCH_DOUBLESIZE &&
                     remote_neighbor && read_parallel_patches)
            {
                /* The coarse grid is now the remote patch;  swap contexts and 
                call same routines above, but with remote patch as the "coarse" 
                grid */
                if (average_from_neighbor)
                {
                    /* Average from local fine grid to remote coarse ghost, since interpolation
                    stencil may need corner coarse grid values */
                    fclaw_patch_average_corner(s->glob,
                                               tdata_fine.this_patch,
                                               tdata_fine.neighbor_patch,
                                               tdata_fine.this_blockno,
                                               tdata_fine.neighbor_blockno,
                                               tdata_fine.is_block_corner,
                                               tdata_fine.icorner,
                                               time_interp,
                                               &tdata_fine);                        

                }
                else if (interpolate_to_neighbor) 
                {
                    /* Interpolate from remote coarse grid patch (coarse grid) to
                       local fine grid patch.  We do not need to average to the 
                       remote patch corners unless corners are used in the 
                       interpolation stencil. */
                       fclaw_patch_interpolate_corner(s->glob,
                                                      tdata_fine.this_patch,
                                                      tdata_fine.neighbor_patch,
                                                      tdata_fine.this_blockno,
                                                      tdata_fine.neighbor_blockno,
                                                      tdata_fine.is_block_corner,
                                                      tdata_fine.icorner,
                                                      time_interp,
                                                      &tdata_fine);

                }
            } /* End of parallel case */
        }  /* End of 'interior_corner' */
    }  /* End of icorner loop */
}

void fclaw_corner_neighbor_indirect(struct fclaw_global* glob,
                                    int minlevel,
                                    int maxlevel,
                                    int time_interp)
{
    const int num_faces = fclaw_domain_num_faces(glob->domain);
    const int num_corners = fclaw_domain_num_corners(glob->domain);

    for(int i = 0; i < glob->domain->num_ghost_patches; i++)
    {
        fclaw_patch_t* this_patch = &glob->domain->ghost_patches[i];
        int level = this_patch->level;
        if (level < minlevel)
        {
            /* We don't need to worry about ghost patches that are at
               coarser levels than we are currently working on */
            continue;
        }

        int this_blockno = fclaw_patch_get_ghost_block(this_patch);

        int intersects_block[num_faces];
        fclaw_block_get_block_boundary(glob, this_patch, intersects_block);

        int intersects_bdry[num_faces];
        for(int iface = 0; iface < num_faces; iface++)
        {
            intersects_bdry[iface] = intersects_block[iface] && glob->domain->blocks[this_blockno].is_boundary[iface];
        }

        /* Transform data needed at multi-block boundaries */
        fclaw_patch_transform_data_t tdata;

        tdata.glob = glob;
        tdata.based = 1;   // cell-centered data in this routine.
        tdata.this_patch = this_patch;
        tdata.this_blockno = this_blockno;
        tdata.this_patchno = i;
        tdata.neighbor_patch = NULL;  // gets filled in below.

        fclaw_patch_transform_init_data(glob,this_patch,
                                        this_blockno,
                                        i,
                                        &tdata);


        fclaw_patch_transform_data_t tdata_fine;

        tdata_fine.glob = glob;
        tdata_fine.neighbor_patch = this_patch;
        tdata_fine.neighbor_blockno = this_blockno;
        tdata_fine.neighbor_patchno = i;
        tdata_fine.based = 1;   // cell-centered data in this routine.

        fclaw_patch_transform_init_data(glob,this_patch,
                                        this_blockno,
                                        i,
                                        &tdata_fine);
        for (int icorner = 0; icorner < num_corners; icorner++)
        {
            int block_corner_count = 0;
            // get corner type and initialize is_block_* values in transform_data and transform_data_finegrid
            // block_i* values will be initialized for transform_data
            // block_i* values will be initialized to -1 for transform_data_finegrid
            int is_interior_in_domain;
            get_corner_type(glob,icorner,
                            intersects_bdry,
                            intersects_block,
                            &is_interior_in_domain,
                            &tdata,
                            &tdata_fine);

            /* Sets block_corner_count to 0 */
            fclaw_patch_set_block_corner_count(glob, this_patch,
                                               icorner,block_corner_count);

            if (is_interior_in_domain)
            {
                /* Is an interior patch corner;  may also be a block corner */

                int is_valid_neighbor;

                int indirect = 1;
                get_corner_neighbor(indirect,
                                    glob,
                                    icorner,
                                    &is_valid_neighbor,
                                    &block_corner_count,
                                    &tdata,
                                    &tdata_fine);

                if (!is_valid_neighbor)
                {
                    /* No corner neighbor.  Either :
                       -- Hanging node
                       -- Cubed sphere
                    */
                    continue;
                }

                if (tdata.neighbor_type == FCLAW_PATCH_SAMESIZE)
                {
                    /* Copy from same size neighbor */
                       fclaw_patch_copy_corner(glob,
                                               tdata.this_patch,
                                               tdata.neighbor_patch,
                                               tdata.this_blockno,
                                               tdata.neighbor_blockno,
                                               tdata.is_block_corner,
                                               tdata.icorner, 
                                               time_interp,
                                               &tdata);
                    ++glob->count_multiproc_corner;
                }
                else if (tdata.neighbor_type == FCLAW_PATCH_HALFSIZE)
                {
                    /* Average from fine grid neighbor */
                       fclaw_patch_average_corner(glob,
                                                  tdata.this_patch,
                                                  tdata.neighbor_patch,
                                                  tdata.this_blockno,
                                                  tdata.neighbor_blockno,
                                                  tdata.is_block_corner,
                                                  tdata.icorner,
                                                  time_interp,
                                                  &tdata);                        
                    ++glob->count_multiproc_corner;
                }
                else if (tdata.neighbor_type == FCLAW_PATCH_DOUBLESIZE)
                {
                    /* Don't do anything; we don't need fine grid ghost cells
                       on ghost patches.  Proof : Consider the corners of the fine
                       patch at either end of the face shared by the coarse and
                       fine patch. Well-balancing assures that at neither of these
                       corners is the fine grid a "coarse grid" to a corner adjacent
                       patch.  So the fine grid will never be needed for interpolation
                       at any grid adjacent to either of these two corners, and so
                       it does not need valid ghost cells along the face shared with the
                       coarse grid. */
                }
            }
        }  /* End of icorner loop */
    }
}