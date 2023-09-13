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

#include <fclaw_edge_neighbors.h>

#include <fclaw_block.h>
#include <fclaw_ghost_fill.h>
#include <fclaw2d_map_query.h>
#include <fclaw_options.h>
#include <fclaw_global.h>
#include <fclaw_physical_bc.h>
#include <fclaw_patch.h>

#include <forestclaw.h>
#include <fclaw_convenience.h>


/* This is used to determine neighbor patch relative level (finer, coarser or samesize)
   This enum is defined both here and in fclaw2d_face_neighbors.cpp.  Is that okay? */
enum
{
    COARSER_GRID = -1,
    SAMESIZE_GRID,
    FINER_GRID
};


static
int get_num_intersections(int intersects[], int faces[])
{
    int num_intersections = 0;
    int i;
    for (i = 0; i < 2; i++)
    {
        if (intersects[faces[i]])
        {
            num_intersections++;
        }
    }
    return num_intersections;
}

static
int find_face(int intersects[], int faces[])
{
    int i;
    for (i = 0; i < 2; i++)
    {
        if (intersects[faces[i]])
        {
            return faces[i];
        }
    }
    return -1;
}

static
void get_edge_type(fclaw_global_t* glob,
                   int iedge,
                   int intersects_bdry[],
                   int intersects_block[],
                   int *interior_edge,
                   int *is_block_edge,
                   int *block_iface)
{
    fclaw_domain_t *domain = glob->domain;

    // p4est has tons of lookup table like this, can be exposed similarly
    int edge_faces[2];
    fclaw_domain_edge_faces(domain, iedge, edge_faces);

    /* Both faces are at a physical boundary */
    int num_phys_faces = get_num_intersections(intersects_bdry,
                                               edge_faces);

    /* Either a edge is at a block boundary (but not a physical boundary),
       or internal to a block.  L-shaped domains are excluded for now
       (i.e. no reentrant corners). */
    *interior_edge = num_phys_faces == 0;

    int num_block_faces = get_num_intersections(intersects_block,
                                                edge_faces);
    /* Both faces are at a block boundary, physical or not */
    *is_block_edge = num_block_faces == domain->dim;

    *block_iface = -1;
    if (num_block_faces == 1)
    {
        *block_iface = find_face(intersects_block, edge_faces);
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
void get_edge_neighbors(fclaw_global_t *glob,
                        int this_block_idx,
                        int this_patch_idx,
                        fclaw_patch_t* this_patch,
                        int icorner,
                        int block_iface,
                        int is_block_edge,
                        int *edge_block_idx,
                        fclaw_patch_t* edge_patches[],
                        int *rcornerno,
                        int **ref_flag_ptr,
                        int *block_edge_count,
                        int ftransform[],
                        fclaw_patch_transform_data_t* ftransform_finegrid)
{
    fclaw_domain_t *domain = glob->domain;
    const int num_face_neighbors = fclaw_domain_num_children(domain)/2;

    for(int i=0; i<2; i++)
    {
        edge_patches[i] = NULL;
    }

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);

#if 0
    /* Note : Pillowsphere case does not return a block corner neighbor */
    int ispillowsphere = 0;
    if(domain->dim == 2)
    {
        //TODO 3d
        fclaw2d_map_pillowsphere(glob);
    }
#endif

    /* See what p4est thinks we have for edges, and consider four cases */
    int rproc_edge[2];
    int edge_patch_idx[2];
    fclaw_patch_relation_t neighbor_type;
    int has_edge_neighbor =
        fclaw_patch_edge_neighbors(domain,
                                   this_block_idx,
                                   this_patch_idx,
                                   icorner,
                                   rproc_edge,
                                   edge_block_idx,
                                   edge_patch_idx,
                                   rcornerno,
                                   &neighbor_type);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);    

    *block_edge_count = 0;  /* Assume we are not at a block edge */
    if (has_edge_neighbor && is_block_edge)
    {
        /* Case 1 : 4 or more patches meet at a block corner.
        This case does NOT include the pillowgrid.   */
        *block_edge_count = 4;  /* assume four for now */

        /* No block corner transforms yet, so we use the 
        interior 'default' transforms. */
        fclaw_patch_transform_blockface_intra (glob, ftransform);
        fclaw_patch_transform_blockface_intra
            (glob, ftransform_finegrid->transform);
    }
    else if (!has_edge_neighbor && !is_block_edge)
    {
        /* Case 2 : 'iedge' is a hanging node */
        /* We do not return valid transformation objects! */
        *ref_flag_ptr = NULL;
        return;
    }
    else if (has_edge_neighbor && !is_block_edge)
    {
        /* Case 3 : 'iedge' is an interior edge, at a block edge,
         or we are on a periodic block.  Need to return a valid
         transform in 'ftransform' */
        /* block_iface is the block number at the of the neighbor? */
        if (block_iface >= 0)
        {
            /* The corner is on a block face (but is not a block edge).
               Compute a transform between blocks. First, get the
               remote face number.  The remote face number encodes the
               orientation, so we have 0 <= rfaceno < 8 */
            int rfaceno;
            int rproc[num_face_neighbors];
            int rpatchno[num_face_neighbors];
            int rblockno;  /* Should equal *corner_block_idx, above. */
            fclaw_patch_face_neighbors(domain,
                                       this_block_idx,
                                       this_patch_idx,
                                       block_iface,
                                       rproc,
                                       &rblockno,
                                       rpatchno,
                                       &rfaceno);

            FCLAW_ASSERT(rblockno == *edge_block_idx);

            /* Get encoding of transforming a neighbor coordinate across a face */
            fclaw_patch_transform_blockface (glob, block_iface, rfaceno, ftransform);

            /* Get transform needed to swap parallel ghost patch with fine
               grid on-proc patch.  This is done so that averaging and
               interpolation routines can be re-used. */
            int iface1 = block_iface;
            int rface1 = rfaceno;
            fclaw_patch_face_swap(domain->dim, &iface1, &rface1);
            fclaw_patch_transform_blockface(glob, iface1, rface1,
                                              ftransform_finegrid->transform);

            ftransform_finegrid->block_iface = iface1;
        }
        else if (this_block_idx == *edge_block_idx)
        {
            /* Both patches are in the same block, so we set the transform to
               a default transform.  This could be the case for periodic boundaries. */
            *block_edge_count = fclaw_domain_num_edges(domain);  /* assume four for now */
            fclaw_patch_transform_blockface_intra (glob, ftransform);
            fclaw_patch_transform_blockface_intra
                (glob, ftransform_finegrid->transform);

        }
        else
        {
            fclaw_global_essentialf("WARNING : this_block_idx = %d\n",this_block_idx);
            fclaw_global_essentialf("WARNING : corner_block_idx = %d\n",*edge_block_idx);
            fclaw_global_essentialf("get_corner_neighbors " \
                                    "(fclaw2d_corner_neighbors.c : " \
                                    "We should not be here\n");
            exit(0);
        }
    }
    else if (!has_edge_neighbor && is_block_edge)
    {
        /* Case 4 : Pillow sphere case or cubed sphere  */
#if 0
        if (!ispillowsphere)
        {
            *block_edge_count = 3;
            /* Exactly 3 patches meet at a corner, e.g. the cubed sphere.
               In this case, 'this_patch' has no corner-adjacent only
               neighbors, and so there is nothing to do. */
            *ref_flag_ptr = NULL;
            *edge_patch = NULL;
            return;
        }
        else
        {
            *block_edge_count = 2;
            has_edge_neighbor = 1;
            int rpatchno[refine_factor];
            int rproc[refine_factor];
            int rfaceno;

            /* Use only faces 0 or 1 to get block data. */
            int iface = icorner % 2;
            neighbor_type =
                fclaw_patch_face_neighbors(domain,
                                             this_block_idx,
                                             this_patch_idx,
                                             iface,
                                             rproc,
                                             edge_block_idx,
                                             rpatchno,
                                             &rfaceno);

            int igrid;
            if (neighbor_type == FCLAW_PATCH_HALFSIZE)
            {
                /* igrid = 0 at corners 0,1 and (R-1) at corners 2,3,
                   where R = refinement factor */
                igrid = (icorner/2)*(refine_factor - 1);
            }
            else
            {
                /* Same size or double size grids. */
                igrid = 0;
            }

            *rcornerno = icorner;    /* This wasn't being set! */
            edge_patch_idx = rpatchno[igrid];
            rproc_edge = rproc[igrid];
        }
#endif
    }

    /* ---------------------------------------------------------------------
       We have a valid neighbor and possibly a transform. We just now need
       to get a pointer to the neighbor patch (which may be a parallel patch)
       and the relative level (-1,0,1).
       --------------------------------------------------------------------- */

    

    int num_neighbors;
    if (neighbor_type == FCLAW_PATCH_HALFSIZE)
    {
        num_neighbors = 2;
        **ref_flag_ptr = 1;
    }
    else if (neighbor_type == FCLAW_PATCH_SAMESIZE)
    {
        num_neighbors = 1;
        **ref_flag_ptr = 0;
    }
    else /* FCLAW2D_PATCH_DOUBLESIZE */
    {
        num_neighbors = 1;
        **ref_flag_ptr = -1;
    }

    for(int i=0; i < num_neighbors; i++)
    {
        if (domain->mpirank != rproc_edge[i])
        {
            edge_patches[i] = &domain->ghost_patches[edge_patch_idx[i]];
        }
        else
        {
            fclaw_block_t *neighbor_block = &domain->blocks[*edge_block_idx];
            edge_patches[i] = &neighbor_block->patches[edge_patch_idx[i]];
        }
    }
}



void cb_edge_fill(fclaw_domain_t *domain,
                  fclaw_patch_t *this_patch,
                  int this_block_idx,
                  int this_patch_idx,
                  void *user)
{
    const int num_faces = fclaw_domain_num_faces(domain);
    const int num_edges = fclaw_domain_num_edges(domain);

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
    int intersects_block[num_faces];

    fclaw_physical_get_bc(s->glob,this_block_idx,this_patch_idx,
                            intersects_bdry);

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


    for (int iedge = 0; iedge < num_edges; iedge++)
    {
        transform_data.iedge = iedge;

        int is_block_edge;
        int is_interior_edge;
        int block_edge_count = 0;
        get_edge_type(s->glob,iedge,
                        intersects_bdry,
                        intersects_block,
                        &is_interior_edge,
                        &is_block_edge,
                        &transform_data.block_iface);

        transform_data_finegrid.block_iface = -1;

        /* Sets block_corner_count to 0 */
        // TODO is this necessary?
        //fclaw_patch_set_block_edge_count(s->glob, this_patch,
        //                                 iedge,block_edge_count);

        if (is_interior_edge)
        {
            /* Is an interior patch edge;  may also be a block edge */

            int edge_block_idx;
            int neighbor_level;
            int *ref_flag_ptr = &neighbor_level;
            fclaw_patch_t* edge_patches[2];
            int redgeno;

            int block_iface = transform_data.block_iface;
            transform_data.is_block_edge = is_block_edge;
            edge_block_idx = -1;
            get_edge_neighbors(s->glob,
                              this_block_idx,
                              this_patch_idx,
                              this_patch,
                              iedge,
                              block_iface,
                              is_block_edge,
                              &edge_block_idx,
                              edge_patches,
                              &redgeno,
                              &ref_flag_ptr,
                              &block_edge_count,
                              transform_data.transform,
                              &transform_data_finegrid);

            /* This sets value in block_corner_count_array */
            //TODO is this necessary?
            //fclaw_patch_set_block_corner_count(s->glob, this_patch,
            //                                     icorner,block_corner_count);

            transform_data.is_block_edge = is_block_edge;

            /* Needed for switching the context */
            transform_data_finegrid.is_block_edge = is_block_edge;
            transform_data_finegrid.iedge = redgeno;
            transform_data_finegrid.this_patch = edge_patches[0];
            transform_data_finegrid.neighbor_patch = this_patch;


            if (ref_flag_ptr == NULL)
            {
                /* No corner neighbor.  Either :
                   -- Hanging node
                   -- Cubed sphere
                */
                continue;
            }

            int remote_neighbor = fclaw_patch_is_ghost(edge_patches[0]);
            if (is_coarse && ((read_parallel_patches && remote_neighbor) || !remote_neighbor))
            {
                transform_data.neighbor_patch = edge_patches[0];
                if (neighbor_level == FINER_GRID)
                {
                    if (interpolate_to_neighbor && !remote_neighbor)
                    {
                        for(int i=0; i < 2; i++)
                        {
                            fclaw_patch_t* coarse_patch = this_patch;
                            fclaw_patch_t* fine_patch = edge_patches[i];
                            transform_data.neighbor_patch = fine_patch;
                            /* No need to interpolate to remote ghost patches. */
                            fclaw_patch_interpolate_edge(s->glob,
                                                         coarse_patch,
                                                         fine_patch,
                                                         iedge,time_interp,
                                                         &transform_data);
                        }
                    }
                    else if (average_from_neighbor)
                    {
                        for(int i=0; i < 2; i++)
                        {
                            /* Average even if neighbor is a remote neighbor */
                            fclaw_patch_t* coarse_patch = this_patch;
                            fclaw_patch_t* fine_patch = edge_patches[i];
                            transform_data.neighbor_patch = fine_patch;
                            fclaw_patch_average_edge(s->glob,
                                                     coarse_patch,
                                                     fine_patch,
                                                     iedge,time_interp,
                                                     &transform_data);                        
                        }
                    }
                }
                else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                {
                    fclaw_patch_copy_edge(s->glob,
                                          this_patch,
                                          edge_patches[0],
                                          this_block_idx,
                                          edge_block_idx,
                                          iedge, time_interp,
                                          &transform_data);
                }

            }  /* End of non-parallel patch case */
            else if (is_fine && neighbor_level == COARSER_GRID &&
                     remote_neighbor && read_parallel_patches)
            {
                /* The coarse grid is now the remote patch;  swap contexts and 
                call same routines above, but with remote patch as the "coarse" 
                grid */
                
                fclaw_patch_t* coarse_patch = edge_patches[0];
                fclaw_patch_t* fine_patch = this_patch;
                int coarse_blockno = edge_block_idx;
                int fine_blockno = this_patch_idx;

                if (interpolate_to_neighbor)
                {
                    /* Interpolate from remote coarse grid patch (coarse grid) to
                       local fine grid patch.  We do not need to average to the 
                       remote patch corners unless corners are used in the 
                       interpolation stencil. */
                    int coarse_iedge = transform_data_finegrid.iedge;

                    fclaw_patch_interpolate_edge(s->glob,
                                                 coarse_patch,
                                                 fine_patch,
                                                 fine_blockno,
                                                 time_interp,
                                                 &transform_data_finegrid);

                }
            } /* End of parallel case */
        }  /* End of 'interior_edge' */
    }  /* End of iedge loop */
}
