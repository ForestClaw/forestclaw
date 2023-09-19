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
                   int *is_interior_in_domain,
                   fclaw_patch_transform_data_t* tdata,
                   fclaw_patch_transform_data_t* tdata_fine)
{
    fclaw_domain_t *domain = glob->domain;

    tdata->icorner = -1;
    tdata->iedge = iedge;
    tdata->iface = -1;

    tdata_fine->icorner = -1;
    tdata_fine->iedge = -1;
    tdata_fine->iface = -1;
    
    //can't be corner
    tdata->is_block_corner = 0;
    tdata_fine->is_block_corner = 0;

    // p4est has tons of lookup table like this, can be exposed similarly
    int edge_faces[2];
    fclaw_domain_edge_faces(domain, iedge, edge_faces);

    /* Both faces are at a physical boundary */
    int num_phys_faces = get_num_intersections(intersects_bdry,
                                               edge_faces);

    /* Either a edge is at a block boundary (but not a physical boundary),
       or internal to a block.  L-shaped domains are excluded for now
       (i.e. no reentrant corners). */
    *is_interior_in_domain = num_phys_faces == 0;

    int num_block_faces = get_num_intersections(intersects_block,
                                                edge_faces);
    /* Both faces are at a block boundary, physical or not */
    tdata->is_block_edge = num_block_faces == domain->refine_dim;
    if(tdata->is_block_edge)
    {
        tdata->block_iedge = iedge;
    }
    else
    {
        tdata->block_iedge = -1;
    }

    tdata->is_block_face = num_block_faces == 1;
    tdata->block_iface = -1;
    if (tdata->is_block_face)
    {
        tdata->block_iface = find_face(intersects_block, edge_faces);
    }

    //set fine tdata
    tdata_fine->is_block_edge = tdata->is_block_edge;
    tdata_fine->is_block_face = tdata->is_block_face;
    tdata_fine->block_iedge = -1;
    tdata_fine->block_iface = -1;
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

typedef struct get_edge_neighbors_return
{
    int is_valid_neighbor;
    fclaw_patch_t* patches[2];
    int patchnos[2];
    int block_edge_count;
} edge_neighbors_t;

static
edge_neighbors_t
get_edge_neighbors(fclaw_global_t *glob,
                   fclaw_patch_transform_data_t* tdata,
                   fclaw_patch_transform_data_t* tdata_nbr)
{
    edge_neighbors_t retval;
    /* assume it is a valid neighbor for now */
    retval.is_valid_neighbor = 1;
    fclaw_domain_t *domain = glob->domain;
    const int num_face_neighbors = fclaw_domain_num_children(domain)/2;

    for(int i=0; i<2; i++)
    {
        retval.patches[i] = NULL;
    }

    fclaw_timer_start (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);

    /* See what p4est thinks we have for edges, and consider four cases */
    int rproc_edge[2];
    int has_edge_neighbor =
        fclaw_patch_edge_neighbors(domain,
                                   tdata->this_blockno,
                                   tdata->this_patchno,
                                   tdata->iedge,
                                   rproc_edge,
                                   &tdata->neighbor_blockno,
                                   retval.patchnos,
                                   &tdata_nbr->iedge,
                                   &tdata->neighbor_type);

    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);    

    retval.block_edge_count = 0;  /* Assume we are not at a block edge */
    if (has_edge_neighbor && tdata->is_block_edge)
    {
        /* Case 1 : 4 or more patches meet at a block corner.
        This case does NOT include the pillowgrid.   */
        retval.block_edge_count = 4;  /* assume four for now */

        /* No block corner transforms yet, so we use the 
        interior 'default' transforms. */
        fclaw_patch_transform_blockface_intra (glob, tdata->transform);
        fclaw_patch_transform_blockface_intra
            (glob, tdata_nbr->transform);
    }
    else if (has_edge_neighbor && tdata->is_block_face)
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
                                   tdata->this_blockno,
                                   tdata->this_patchno,
                                   tdata->block_iface,
                                   rproc,
                                   &rblockno,
                                   rpatchno,
                                   &rfaceno);

        FCLAW_ASSERT(rblockno == tdata->neighbor_blockno);

        /* Get encoding of transforming a neighbor coordinate across a face */
        fclaw_patch_transform_blockface (glob, tdata->block_iface, rfaceno, tdata->transform);

        /* Get transform needed to swap parallel ghost patch with fine
           grid on-proc patch.  This is done so that averaging and
           interpolation routines can be re-used. */
        int iface1 = tdata->block_iface;
        int rface1 = rfaceno;
        fclaw_patch_face_swap(domain->refine_dim, &iface1, &rface1);
        fclaw_patch_transform_blockface(glob, iface1, rface1,
                                          tdata_nbr->transform);

        tdata_nbr->block_iface = iface1;
    }
    else if (has_edge_neighbor)
    {
        /* Case 3 : 'iedge' is an interior edge */
        if (tdata->this_blockno == tdata->neighbor_blockno)
        {
            /* Both patches are in the same block, so we set the transform to
               a default transform.  This could be the case for periodic boundaries. */
            retval.block_edge_count = 4;
            fclaw_patch_transform_blockface_intra(glob, tdata->transform);
            fclaw_patch_transform_blockface_intra(glob, tdata_nbr->transform);

        }
        else
        {
            fclaw_errorf("WARNING : this_block_idx = %d\n",tdata->this_blockno);
            fclaw_errorf("WARNING : neighbor_block_idx = %d\n",tdata->neighbor_blockno);
            fclaw_abortf("get_edge_neighbors " \
                            "(fclaw_edge_neighbors.c : " \
                            "We should not be here\n");
        }
    }
    else if (!has_edge_neighbor)
    {
        /* Case 2 : 'iedge' is a hanging node */
        /* We do not return valid transformation objects! */
        retval.is_valid_neighbor = 0;
        return retval;
    }

    /* ---------------------------------------------------------------------
       We have a valid neighbor and possibly a transform. We just now need
       to get a pointer to the neighbor patch (which may be a parallel patch)
       and the relative level (-1,0,1).
       --------------------------------------------------------------------- */

    
    int num_neighbors;
    if (tdata->neighbor_type == FCLAW_PATCH_HALFSIZE)
    {
        tdata_nbr->neighbor_type = FCLAW_PATCH_DOUBLESIZE;
        num_neighbors = 2;
    }
    else if (tdata->neighbor_type == FCLAW_PATCH_SAMESIZE)
    {
        tdata_nbr->neighbor_type = FCLAW_PATCH_SAMESIZE;
        num_neighbors = 1;
    }
    else /* FCLAW2D_PATCH_DOUBLESIZE */
    {
        tdata_nbr->neighbor_type = FCLAW_PATCH_HALFSIZE;
        num_neighbors = 1;
    }

    for(int i=0; i < num_neighbors; i++)
    {
        if (domain->mpirank != rproc_edge[i])
        {
            retval.patches[i] = &domain->ghost_patches[retval.patchnos[i]];
        }
        else
        {
            fclaw_block_t *neighbor_block = &domain->blocks[tdata->neighbor_blockno];
            retval.patches[i] = &neighbor_block->patches[retval.patchnos[i]];
        }
    }

    //set tdata_fine
    tdata_nbr->this_blockno = tdata->neighbor_blockno;

    return retval;
}



void cb_edge_fill(fclaw_domain_t *domain,
                  fclaw_patch_t *this_patch,
                  int this_blockno,
                  int this_patchno,
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

    fclaw_physical_get_bc(s->glob,this_blockno,this_patchno,
                            intersects_bdry);

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


    fclaw_patch_transform_data_t tdata_nbr;

    tdata_nbr.glob = s->glob;
    tdata_nbr.based = 1;   // cell-centered data in this routine.
    tdata_nbr.neighbor_patch = this_patch;
    tdata_nbr.neighbor_blockno = this_blockno;
    tdata_nbr.neighbor_patchno = this_patchno;

    fclaw_patch_transform_init_data(s->glob,this_patch,
                                      this_blockno,
                                      this_patchno,
                                      &tdata_nbr);


    for (int iedge = 0; iedge < num_edges; iedge++)
    {
        int is_interior_in_domain;
        get_edge_type(s->glob,iedge,
                        intersects_bdry,
                        intersects_block,
                        &is_interior_in_domain,
                        &tdata,
                        &tdata_nbr);

        /* Sets block_corner_count to 0 */
        // TODO is this necessary?
        //fclaw_patch_set_block_edge_count(s->glob, this_patch,
        //                                 iedge,block_edge_count);

        if (is_interior_in_domain)
        {
            /* Is an interior patch edge;  may also be a block edge */

            edge_neighbors_t edge_neighbors =
                get_edge_neighbors(s->glob,
                                   &tdata,
                                   &tdata_nbr);

            /* This sets value in block_corner_count_array */
            //TODO is this necessary?
            //fclaw_patch_set_block_corner_count(s->glob, this_patch,
            //                                     icorner,block_corner_count);

            if (!edge_neighbors.is_valid_neighbor)
            {
                /* No corner neighbor.  Either :
                   -- Hanging node
                   -- Cubed sphere
                */
                continue;
            }
            FCLAW_ASSERT(tdata.iedge == iedge);
            FCLAW_ASSERT(tdata_nbr.iedge == (iedge^3));

            /* Parallel distribution keeps siblings on same processor */
            int remote_neighbor = fclaw_patch_is_ghost(edge_neighbors.patches[0]);
            if (is_coarse)
            {
                if (tdata.neighbor_type == FCLAW_PATCH_HALFSIZE)
                {
                    for (int igrid = 0; igrid < 2; igrid++)
                    {
                        tdata.neighbor_patch = edge_neighbors.patches[igrid];
                        tdata.neighbor_patchno = edge_neighbors.patchnos[igrid];

                        remote_neighbor = fclaw_patch_is_ghost(tdata.neighbor_patch);
                        int valid_remote = read_parallel_patches && remote_neighbor;
                        int local_neighbor = !remote_neighbor;
                        if (!(local_neighbor || valid_remote))
                        {
                            continue;
                        }

                        if (interpolate_to_neighbor && !remote_neighbor)
                        {
                            /* interpolate to igrid */
                            fclaw_patch_interpolate_edge(s->glob,
                                                         tdata.this_patch,
                                                         tdata.neighbor_patch,
                                                         tdata.iedge,
                                                         time_interp,
                                                         &tdata);
                        }
                        else if (average_from_neighbor)
                        {
                            /* average from igrid */
                            fclaw_patch_average_edge(s->glob,
                                                     tdata.this_patch,
                                                     tdata.neighbor_patch,
                                                     tdata.iedge,
                                                     time_interp,
                                                     &tdata);
                        }
                    }
                }
                else if (tdata.neighbor_type == FCLAW_PATCH_SAMESIZE && copy_from_neighbor)
                {
                    /* Copy to same size patch */
                    tdata.neighbor_patch = edge_neighbors.patches[0];
                    tdata.neighbor_patchno = edge_neighbors.patchnos[0];

                    fclaw_patch_copy_edge(s->glob,
                                          tdata.this_patch,
                                          tdata.neighbor_patch,
                                          tdata.this_blockno,
                                          tdata.neighbor_blockno,
                                          tdata.iedge,
                                          time_interp,
                                          &tdata);

                    /* We also need to copy _to_ the remote neighbor; switch contexts, but
                       use ClawPatches that are only in scope here, to avoid
                       conflicts with above uses of the same variables. This is needed
                       in case we want to interpolate to adjacent corners on fine grids.*/
                    if (remote_neighbor)
                    {
                        tdata_nbr.this_patch = edge_neighbors.patches[0];
                        tdata_nbr.this_patchno = edge_neighbors.patchnos[0];
                        fclaw_patch_copy_edge(s->glob,
                                              tdata_nbr.this_patch,
                                              tdata_nbr.neighbor_patch,
                                              tdata_nbr.this_blockno,
                                              tdata_nbr.neighbor_blockno,
                                              tdata_nbr.iedge, 
                                              time_interp,
                                              &tdata_nbr);                            
                    }
                }
            }
            else if (is_fine && tdata.neighbor_type == FCLAW_PATCH_DOUBLESIZE && remote_neighbor
                     && read_parallel_patches)
            {
                /* Swap 'this_patch' (fine grid) and the neighbor patch 
                (a coarse grid) */
                tdata_nbr.this_patch = edge_neighbors.patches[0];
                tdata_nbr.this_patchno = edge_neighbors.patchnos[0];

                if (average_from_neighbor)
                {
                        /* Average from 'this' grid (fine grid) to remote grid 
                    (coarse grid) */
                    fclaw_patch_average_edge(s->glob,
                                             tdata_nbr.this_patch,
                                             tdata_nbr.neighbor_patch,
                                             tdata_nbr.iedge,
                                             time_interp,
                                             &tdata_nbr);
                }
                else if (interpolate_to_neighbor)
                {
                    /* Interpolate from remote neighbor to 'this' patch (the finer grid */
                    fclaw_patch_interpolate_edge(s->glob,
                                                 tdata_nbr.this_patch,
                                                 tdata_nbr.neighbor_patch,
                                                 tdata_nbr.iedge,
                                                 time_interp,
                                                 &tdata_nbr);
                }
            }
        }  /* End of 'interior_edge' */
    }  /* End of iedge loop */
}
