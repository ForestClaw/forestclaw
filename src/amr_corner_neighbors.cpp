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

#include "amr_includes.H"

void
get_corner_type(fclaw2d_domain_t* domain,
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

    // Both faces are at a physical boundary
    fclaw_bool is_phys_corner =
        intersects_bdry[corner_faces[0]] && intersects_bdry[corner_faces[1]];

    // Corner lies in interior of physical boundary edge.
    fclaw_bool corner_on_phys_face = !is_phys_corner &&
             (intersects_bdry[corner_faces[0]] || intersects_bdry[corner_faces[1]]);

    /* Either a corner is at a block boundary (but not a physical boundary),
       or internal to a block.  L-shaped domains are excluded for now
       (i.e. no reentrant corners). */
    *interior_corner = !corner_on_phys_face && !is_phys_corner;

    // Both faces are at a block boundary
    *is_block_corner =
        intersects_block[corner_faces[0]] && intersects_block[corner_faces[1]];

    *block_iface = -1;
    if (!*is_block_corner)
    {
        if (intersects_block[corner_faces[0]])
        {
            // Corner is on a block face.
            *block_iface = corner_faces[0];
        }
        else if (intersects_block[corner_faces[1]])
        {
            // Corner is on a block face.
            *block_iface = corner_faces[1];
        }
    }
}


/* -----------------------------------------------
   Four cases to consider.   The 'has_corner_neighbor'
   value is returned from p4est.

   Case No. | has_corner_neighbor  |  is_block_corner
   ----------------------------------------------------
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
   ----------------------------------------------- */


void get_corner_neighbor(fclaw2d_domain_t *domain,
                         int this_block_idx,
                         int this_patch_idx,
                         int icorner,
                         int iface,
                         fclaw_bool is_block_corner,
                         int *corner_block_idx,
                         fclaw2d_patch_t** ghost_patch,
                         int **ref_flag_ptr,
                         int ftransform[])
{
    /* See what p4est thinks we have for corners, and consider four cases */
    int rproc_corner;
    int corner_patch_idx;
    fclaw2d_patch_relation_t neighbor_type;

    fclaw_bool has_corner_neighbor =
        fclaw2d_patch_corner_neighbors(domain, this_block_idx, this_patch_idx,
                                       icorner, &rproc_corner, corner_block_idx,
                                       &corner_patch_idx, &neighbor_type);

    if (has_corner_neighbor && is_block_corner)
    {
        /* 4 or more patches meet at a block corner */
        printf("get_corner_neighbors (amr_neighbors.f) : "             \
               "Four or more patches meet at a block corner; "         \
               "(not implemented yet).\n");
        exit(0);
    }
    else if (!has_corner_neighbor && !is_block_corner)
    {
        /* 'icorner' is a hanging node */
        *ref_flag_ptr = NULL;
        *ghost_patch = NULL;
        return;
    }
    else if (has_corner_neighbor && !is_block_corner)
    {
        /* 'icorner' is an interior corner or a a block edge */
        /* Need to return a valid transform in 'ftransform' */
        if (this_block_idx == *corner_block_idx)
        {
            /* Both patches are in the same block, so we set the transform to
               a default transform */
            ftransform[8] = 4;
        }
        else
        {
            /* The corner is on a block edge (but is not a block corner).
               In this case, we need to compute a transform between blocks.
               To to this, we need the remote face number */
            int rfaceno;
            int rproc[p4est_refineFactor];
            int rpatchno;
            int rblockno;  /* Should equal *corner_block_idx, above. */
            fclaw2d_patch_face_neighbors(domain,
                                         this_block_idx,
                                         this_patch_idx,
                                         iface,
                                         rproc,
                                         &rblockno,
                                         &rpatchno,
                                         &rfaceno);

            /* Get encoding of transforming a neighbor coordinate across a face */
            fclaw2d_patch_face_transformation (iface, rfaceno, ftransform);
        }
    }
    else if (!has_corner_neighbor && is_block_corner)
    {
        if (!ispillowsphere_())
        {
            /* Exactly 3 patches meet at a corner, e.g. the cubed sphere.
               In this case, 'this_patch' has no corner-adjacent only
               neighbors, and so there is nothing to do. */
            *ref_flag_ptr = NULL;
            *ghost_patch = NULL;
            return;
        }
        else
        {
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
       to get a pointer to the ghost patch (which may be a parallel patch)
       and the relative level (-1,0,1).
       --------------------------------------------------------------------- */

    if (domain->mpirank != rproc_corner)
    {
        *ghost_patch = &domain->ghost_patches[corner_patch_idx];
    }
    else
    {
        fclaw2d_block_t *neighbor_block = &domain->blocks[*corner_block_idx];
        *ghost_patch = &neighbor_block->patches[corner_patch_idx];
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
    fclaw_bool is_coarse = filltype->is_coarse;
    fclaw_bool is_fine = filltype->is_fine;

    fclaw_bool intersects_bdry[NumFaces];
    fclaw_bool intersects_block[NumFaces];
    fclaw_bool is_block_corner;
    fclaw_bool is_interior_corner;

    get_phys_boundary(domain,this_block_idx,this_patch_idx,
                      intersects_bdry);

    get_block_boundary(domain,this_block_idx,this_patch_idx,
                       intersects_block);

    /* Transform data needed at multi-block boundaries */
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;   // cell-centered data in this routine.
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;  // gets filled in below.

    int refratio = gparms->refratio;

    for (int icorner = 0; icorner < NumCorners; icorner++)
    {
        get_corner_type(domain,icorner,
                        intersects_bdry,
                        intersects_block,
                        &is_interior_corner,
                        &is_block_corner,
                        &transform_data.iface);

        if (is_interior_corner)  /* Interior to the domain, not necessarily to a block */
        {
            int corner_block_idx;
            int relative_refinement_level;
            int *ref_flag_ptr = &relative_refinement_level;
            fclaw2d_patch_t *ghost_patch;

            transform_data.icorner = icorner;
            get_corner_neighbor(domain,
                                this_block_idx,
                                this_patch_idx,
                                icorner,
                                transform_data.iface,
                                is_block_corner,
                                &corner_block_idx,
                                &ghost_patch,
                                &ref_flag_ptr,
                                transform_data.transform);

            if (ref_flag_ptr == NULL)
            {
                /* no corner neighbor; relative_refinement_level is not set
                   This can happen in the cubed sphere case, or if icorner is
                   a hanging node */
                continue;
            }

            if (filltype->copy)
            {
                if (relative_refinement_level == 0)
                {
                    ClawPatch *this_cp = get_clawpatch(this_patch);
                    fclaw2d_patch_t* corner_patch = ghost_patch;
                    ClawPatch *corner_cp = get_clawpatch(corner_patch);
                    transform_data.neighbor_patch = corner_patch;
                    if (!is_block_corner)
                    {
                        this_cp->exchange_corner_ghost(icorner,corner_cp,
                                                       &transform_data);
                    }
                    else
                    {
                        if (ispillowsphere_())
                        {
                            this_cp->mb_exchange_block_corner_ghost(icorner,corner_cp);
                        }
                        else
                        {
                            /* Handle 4 and 5 corner block cases here;  nothing to do for
                               cubed sphere case. */
                        }
                    }
                }
            }
            else if (filltype->average)
            {
                if (relative_refinement_level == 1 && is_coarse)
                {
                    /* Corner neighbor at a finer level, and so we need to average
                       that corner onto the coarser corner ghost cells */
                    ClawPatch *this_cp = get_clawpatch(this_patch);
                    ClawPatch *corner_cp = get_clawpatch(ghost_patch);
                    transform_data.neighbor_patch = ghost_patch;

                    if (!is_block_corner)
                    {
                        this_cp->average_corner_ghost(icorner,refratio,corner_cp,
                                                      time_interp, &transform_data);
                    }
                    else
                    {
                        if (ispillowsphere_())
                        {
                            this_cp->mb_average_block_corner_ghost(icorner,refratio,
                                                                   corner_cp,time_interp);
                        }
                        else
                        {
                            /* Handle 4 and 5 corner block cases here */
                        }
                    }
                }
                else if (relative_refinement_level == -1 && is_fine)
                {
                    /* Neighbor is a parallel patch;  swap 'this' and 'neighbor' */
                }
            }
            else if (filltype->interpolate)
            {
                if (relative_refinement_level == 1 && is_coarse)
                {
                    ClawPatch *this_cp = get_clawpatch(this_patch);
                    ClawPatch *corner_cp = get_clawpatch(ghost_patch);
                    transform_data.neighbor_patch = ghost_patch;
                    if (!is_block_corner)
                    {
                        this_cp->interpolate_corner_ghost(icorner,refratio,corner_cp,
                                                          time_interp,&transform_data);
                    }
                    else
                    {
                        if (ispillowsphere_())
                        {
                            this_cp->mb_interpolate_block_corner_ghost(icorner,refratio,
                                                                       corner_cp,time_interp);
                        }
                        else
                        {
                            /* Handle 4,5 patch corners; nothing to be done for cubed sphere */
                        }
                    }
                }
                else if (relative_refinement_level == -1 && is_fine)
                {
                    /* Neighbor is a parallel patch;  swap 'this' and 'neighbor' */
                }
            }
        }
    }
}
