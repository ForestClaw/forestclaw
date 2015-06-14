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

#include <fclaw2d_global.h>

#include "fclaw2d_ghost_fill.h"
#include "ClawPatch.H"
#include "amr_includes.H"
#include "fclaw2d_map_query.h"

/* This is used to determine neighbor patch relative level (finer, coarser or samesize)
   This enum is defined both here and in fclaw2d_face_neighbors.cpp.  Is that okay? */
enum
{
    COARSER_GRID = -1,
    SAMESIZE_GRID,
    FINER_GRID
};


static
void get_corner_type(fclaw2d_domain_t* domain,
                     int icorner,
                     fclaw_bool intersects_bdry[],
                     fclaw_bool intersects_block[],
                     fclaw_bool *interior_corner,
                     fclaw_bool *is_block_corner,
                     int *block_iface)
{
    // p4est has tons of lookup table like this, can be exposed similarly
    int corner_faces[SpaceDim];
    fclaw2d_domain_corner_faces(domain, icorner, corner_faces);

    /* Both faces are at a physical boundary */
    fclaw_bool is_phys_corner =
        intersects_bdry[corner_faces[0]] && intersects_bdry[corner_faces[1]];

    /* Corner lies in interior of physical boundary edge. */
    fclaw_bool corner_on_phys_face = !is_phys_corner &&
             (intersects_bdry[corner_faces[0]] || intersects_bdry[corner_faces[1]]);

    /* Either a corner is at a block boundary (but not a physical boundary),
       or internal to a block.  L-shaped domains are excluded for now
       (i.e. no reentrant corners). */
    *interior_corner = !corner_on_phys_face && !is_phys_corner;

    /* Both faces are at a block boundary, physical or not */
    *is_block_corner =
        intersects_block[corner_faces[0]] && intersects_block[corner_faces[1]];

    *block_iface = -1;
    if (!*is_block_corner)
    {
        /* At most one of these is true if corner is not a block corner */
        if (intersects_block[corner_faces[0]])
        {
            /* Corner is on a block face. */
            *block_iface = corner_faces[0];
        }
        else if (intersects_block[corner_faces[1]])
        {
            /* Corner is on a block face. */
            *block_iface = corner_faces[1];
        }
    }
}


/* --------------------------------------------------------
   Four cases to consider.   The 'has_corner_neighbor'
   value is returned from p4est.

   Case No. | has_corner_neighbor  |  is_block_corner
   --------------------------------------------------------
      1     |       T              |        T
      2     |       F              |        F
      3     |       T              |        F
      4     |       F              |        T

    Case 1 : In this case, 4 or more patches meet at a block
             corner. (not implemented yet).
    Case 2 : Corner is a hanging node and has no valid
             adjacent corner.
    Case 3 : Corner is either interior to a block, or on a
             block edge.  In each sub-case, the transform is
             well-defined.
    Case 4 : Either 3 patches meet at a corner, in which
             case we don't have a valid corner, or we are
             on a pillow grid, in which case we have a valid
             corner, but one which we nonetheless treat
             as a special case.
   ------------------------------------------------------ */

static
void get_corner_neighbor(fclaw2d_domain_t *domain,
                         int this_block_idx,
                         int this_patch_idx,
                         fclaw2d_patch_t* this_patch,
                         int icorner,
                         int block_iface,
                         fclaw_bool is_block_corner,
                         int *corner_block_idx,
                         fclaw2d_patch_t** corner_patch,
                         int *rcornerno,
                         int **ref_flag_ptr,
                         int *block_corner_count,
                         int ftransform[])
{
    /* See what p4est thinks we have for corners, and consider four cases */
    int rproc_corner;
    int corner_patch_idx;
    fclaw2d_patch_relation_t neighbor_type;

    fclaw2d_map_context_t *cont = get_map_context(domain);
    fclaw_bool ispillowsphere = FCLAW2D_MAP_IS_PILLOWSPHERE(&cont) != 0; //

    fclaw_bool has_corner_neighbor =
        fclaw2d_patch_corner_neighbors(domain, this_block_idx, this_patch_idx,
                                       icorner, &rproc_corner, corner_block_idx,
                                       &corner_patch_idx, rcornerno,
                                       &neighbor_type);

    if (domain->mpirank == 0 && this_patch_idx == 11 && this_block_idx == 0)
    {
        fclaw_global_essentialf("corner_block_idx = %d\n",*corner_block_idx);
        fclaw_global_essentialf("this_block_idx = %d\n",this_block_idx);
        fclaw_global_essentialf("local_num_patches : %d\n",domain->local_num_patches);
        fclaw_global_essentialf("global_num_patches : %lld\n",domain->global_num_patches);
        fclaw_global_essentialf("ghost_num_patches : %d\n",domain->num_ghost_patches);
    }

    *block_corner_count = 0;  /* Assume we are not at a block corner */
    if (has_corner_neighbor && is_block_corner)
    {
        /* 4 or more patches meet at a block corner */
        *block_corner_count = 4;  /* assume four for now */
        /* We don't have a block corner transformation, so I am going to
           treat this as if it were an interior corner */
        ftransform[8] = 4;
    }
    else if (!has_corner_neighbor && !is_block_corner)
    {
        /* 'icorner' is a hanging node */
        *ref_flag_ptr = NULL;
        *corner_patch = NULL;
        return;
    }
    else if (has_corner_neighbor && !is_block_corner)
    {
        /* 'icorner' is an interior corner, at a block edge,
         or we are on a periodic block.  Need to return a valid
         transform in 'ftransform' */
        if (block_iface >= 0)
        {
            /* The corner is on a block edge (but is not a block corner).
               Compute a transform between blocks. First, get the
               remote face number */
            int rfaceno;
            int rproc[p4est_refineFactor];
            int rpatchno;
            int rblockno;  /* Should equal *corner_block_idx, above. */
            fclaw2d_patch_face_neighbors(domain,
                                         this_block_idx,
                                         this_patch_idx,
                                         block_iface,
                                         rproc,
                                         &rblockno,
                                         &rpatchno,
                                         &rfaceno);

            /* Get encoding of transforming a neighbor coordinate across a face */
            fclaw2d_patch_face_transformation (block_iface, rfaceno, ftransform);
        }
        else if (this_block_idx == *corner_block_idx)
        {
            /* Both patches are in the same block, so we set the transform to
               a default transform */
            ftransform[8] = 4;
        }
        else
        {
            printf("get_corner_neighbors (fclaw2d_corner_neighbors.cpp : We should not be here\n");
            exit(0);
        }
    }
    else if (!has_corner_neighbor && is_block_corner)
    {
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
            /* 2 patches meet at a corner, .e.g. pillow sphere.
               This is handled as a special case */
            has_corner_neighbor = fclaw_true;
            int rpatchno[p4est_refineFactor];
            int rproc[p4est_refineFactor];
            int rfaceno;

            /* Use only faces 0 or 1 to get block data. */
            int iface = icorner % 2;
            neighbor_type =
                fclaw2d_patch_face_neighbors(domain,
                                             this_block_idx,
                                             this_patch_idx,
                                             iface,
                                             rproc,
                                             corner_block_idx,
                                             rpatchno,
                                             &rfaceno);

            int igrid;
            if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
            {
                /* igrid = 0 at corners 0,1 and (R-1) at corners 2,3,
                   where R = refinement factor */
                igrid = (icorner/2)*(p4est_refineFactor - 1);
            }
            else
            {
                /* Same size or double size grids. */
                igrid = 0;
            }

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
        fclaw2d_block_t *neighbor_block = &domain->blocks[*corner_block_idx];
        *corner_patch = &neighbor_block->patches[corner_patch_idx];
    }

    if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
    {
        **ref_flag_ptr = 1;
    }
    else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
    {
        **ref_flag_ptr = 0;
    }
    else /* FCLAW2D_PATCH_DOUBLESIZE */
    {
        **ref_flag_ptr = -1;
    }
}



void cb_corner_fill(fclaw2d_domain_t *domain,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    void *user)
{
    fclaw2d_exchange_info_t *filltype = (fclaw2d_exchange_info_t*) user;
    fclaw_bool time_interp = filltype->time_interp;
    fclaw_bool is_coarse = filltype->grid_type == FCLAW2D_IS_COARSE;
    fclaw_bool is_fine = filltype->grid_type == FCLAW2D_IS_FINE;

    fclaw_bool read_parallel_patches = filltype->read_parallel_patches;

    fclaw_bool copy_from_neighbor = filltype->exchange_type == FCLAW2D_COPY;
    fclaw_bool average_from_neighbor = filltype->exchange_type == FCLAW2D_AVERAGE;
    fclaw_bool interpolate_to_neighbor = filltype->exchange_type == FCLAW2D_INTERPOLATE;

    fclaw_bool intersects_bdry[NumFaces];
    fclaw_bool intersects_block[NumFaces];
    fclaw_bool is_block_corner;
    fclaw_bool is_interior_corner;
    int block_corner_count;

    fclaw2d_map_context_t *cont = get_map_context(domain);
    fclaw_bool ispillowsphere = FCLAW2D_MAP_IS_PILLOWSPHERE(&cont) != 0; //

    fclaw2d_get_physical_bc(domain,this_block_idx,this_patch_idx,
                            intersects_bdry);

    fclaw2d_block_get_block_boundary(domain, this_patch, intersects_block);

    /* Transform data needed at multi-block boundaries */
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.
    transform_data.is_block_corner = 0;

    int refratio = gparms->refratio;

    ClawPatch *this_cp = fclaw2d_clawpatch_get_cp(this_patch);
    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        block_corner_count = 0;
        get_corner_type(domain,icorner,
                        intersects_bdry,
                        intersects_block,
                        &is_interior_corner,
                        &is_block_corner,
                        &transform_data.block_iface);

        if (is_interior_corner)
        {
            /* Interior to the domain, not necessarily to a block */

            int corner_block_idx;
            int neighbor_level;
            int *ref_flag_ptr = &neighbor_level;
            fclaw2d_patch_t *corner_patch;
            int rcornerno;

            transform_data.icorner = icorner;
            corner_block_idx = -1;
            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                this_patch,
                                icorner,
                                transform_data.block_iface,
                                is_block_corner,
                                &corner_block_idx,
                                &corner_patch,
                                &rcornerno,
                                &ref_flag_ptr,
                                &block_corner_count,
                                transform_data.transform);

            this_cp->set_block_corner_count(icorner,block_corner_count);
            transform_data.is_block_corner = is_block_corner;

            if (ref_flag_ptr == NULL)
            {
                /* no corner neighbor; neighbor_level is not set
                   This can happen in the cubed sphere case, or if icorner is
                   a hanging node */
                continue;
            }

            fclaw_bool remote_neighbor = fclaw2d_patch_is_ghost(corner_patch);
            if (is_coarse && ((read_parallel_patches && remote_neighbor) || !remote_neighbor))
            {
                ClawPatch *corner_cp = fclaw2d_clawpatch_get_cp(corner_patch);
                transform_data.neighbor_patch = corner_patch;
                if (!(ispillowsphere && is_block_corner))
                {
                    if (neighbor_level == FINER_GRID)
                    {
                        if (interpolate_to_neighbor && !remote_neighbor)
                        {
                            /* Interpolate 'this_cp' (coarse grid) to 'corner_cp' (fine grid)
                               'icorner' is the coarse grid corner. */
                            this_cp->interpolate_corner_ghost(icorner,refratio,corner_cp,
                                                              time_interp,&transform_data);
                        }
                        else if (average_from_neighbor)
                        {
                            /* average 'corner_cp' (fine grid) to 'this_cp' (coarse grid) */
                            this_cp->average_corner_ghost(icorner,refratio,corner_cp,
                                                          time_interp, &transform_data);
                        }
                    }
                    else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                    {
                        this_cp->exchange_corner_ghost(icorner,corner_cp,
                                                       &transform_data);
                    }
                }
                else /* is_block_corner && ispillowsphere */
                {
                    /* Pillowsphere : The block corners of the pillow sphere have to
                       be handled as a special case */
                    if (neighbor_level == FINER_GRID)
                    {
                        if (interpolate_to_neighbor && !remote_neighbor)
                        {
                            this_cp->mb_interpolate_block_corner_ghost(icorner,refratio,
                                                                       corner_cp,time_interp);
                        }
                        else if (average_from_neighbor)
                        {
                            this_cp->mb_average_block_corner_ghost(icorner,refratio,
                                                                   corner_cp,time_interp);
                        }
                    }
                    else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                    {
                        this_cp->mb_exchange_block_corner_ghost(icorner,corner_cp);
                    }
                }
            }  /* Ende of non-parallel patch case */
            else if (is_fine && neighbor_level == COARSER_GRID &&
                     remote_neighbor && read_parallel_patches)
            {
                /* Swap 'this_patch' and the neighbor patch */
                ClawPatch *coarse_cp = fclaw2d_clawpatch_get_cp(corner_patch);
                ClawPatch *fine_cp = this_cp;

                /* Need to switch corners */
                int coarse_icorner = 3-icorner;

                transform_data.this_patch = corner_patch;
                transform_data.neighbor_patch = this_patch;

		/* "this" grid is now the remote patch; average from "this" to on-proc fine grid */
		if (average_from_neighbor)
                {
		    /* Average from remote patch (fine grid) to 'this_patch' (coarse grid) */
		    coarse_cp->average_corner_ghost(coarse_icorner,refratio,
                                                    fine_cp,time_interp,
                                                    &transform_data);
                }
                else if (interpolate_to_neighbor)
                {
                    /* Interpolate from remote patch (coarse grid) to 'this' patch (fine grid) */
                    coarse_cp->interpolate_corner_ghost(coarse_icorner,refratio,fine_cp,
                                                      time_interp,&transform_data);
                }
            }  /* End of parallel case */
        }  /* End of 'interior_corner' */
    }  /* End of icorner loop */
}
