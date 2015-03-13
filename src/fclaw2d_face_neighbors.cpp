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

/** \file fclaw2d_face_neighbors.cpp
 * Average, coarsen and copy between grids at faces.
 **/

#include <fclaw2d_global.h>

#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_clawpatch.hpp>
#include <ClawPatch.H>

/* This is used to determine neighbor patch relative level (finer, coarser or samesize) */
enum
{
    COARSER_GRID = -1,
    SAMESIZE_GRID,
    FINER_GRID
};

/* This function is a bit overkill, but I put it here so the logic in both
   the corner fill and face fill are the same */
static
void get_face_type(fclaw2d_domain_t *domain,
                   int iface,
                   fclaw_bool intersects_phys_bdry[],
                   fclaw_bool intersects_block[],
                   fclaw_bool *is_block_face,
                   fclaw_bool *is_interior_face)
{
    *is_block_face = intersects_block[iface];
    *is_interior_face = !intersects_phys_bdry[iface];
}


static
void get_face_neighbors(fclaw2d_domain_t *domain,
                        int this_block_idx,
                        int this_patch_idx,
                        int iface,
                        int is_block_face,
                        int *neighbor_block_idx,
                        fclaw2d_patch_t* neighbor_patches[],
                        int **ref_flag_ptr,
                        int **fine_grid_pos_ptr,
                        int **iface_neighbor_ptr,
                        int ftransform[])
{
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;
    int num_neighbors;

    for(int ir = 0; ir < p4est_refineFactor; ir++)
    {
        neighbor_patches[ir] = NULL;
    }

    fclaw2d_patch_relation_t neighbor_type =
        fclaw2d_patch_face_neighbors(domain,
                                     this_block_idx,
                                     this_patch_idx,
                                     iface,
                                     rproc,
                                     &rblockno,
                                     rpatchno,
                                     &rfaceno);

    /* ------------------------------
      neighbor_type is one of :
      FCLAW2D_PATCH_BOUNDARY,    -- physical boundary
      FCLAW2D_PATCH_HALFSIZE,
      FCLAW2D_PATCH_SAMESIZE,
      FCLAW2D_PATCH_DOUBLESIZE
      ------------------------------- */

    if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
    {
        /* This case should be excluded by earlier checkes */
        printf("get_face_neighbors (fclaw2d_face_neighbors.cpp) : No patch found\n");
        exit(0);
    }
    else
    {
        *neighbor_block_idx = is_block_face ? rblockno : -1;
        /* Get encoding of transforming a neighbor coordinate across a face */
        fclaw2d_patch_face_transformation (iface, rfaceno, ftransform);

        if (!is_block_face)
        {
            /* If we are within one patch this is a special case */
            FCLAW_ASSERT (*neighbor_block_idx == -1);
            ftransform[8] = 4;
        }

        if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            **ref_flag_ptr = 0;
            *fine_grid_pos_ptr = NULL;
            num_neighbors = 1;
        }
        else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
        {
            **ref_flag_ptr = -1;
            **fine_grid_pos_ptr = rproc[1];    /* Special storage for fine grid info */
            num_neighbors = 1;
        }
        else if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            /* Patch has two neighbors */
            **ref_flag_ptr = 1; /* patches are at one level finer */
            *fine_grid_pos_ptr = NULL;
            num_neighbors = p4est_refineFactor;
        }
        else
        {
            printf ("Illegal fclaw2d_patch_face_neighbors return value\n");
            exit (1);
        }

        for(int ir = 0; ir < num_neighbors; ir++)
        {
            fclaw2d_patch_t *neighbor;
            if (rproc[ir] == domain->mpirank)
            {
                /* neighbor patch is local */
                fclaw2d_block_t *neighbor_block = &domain->blocks[rblockno];
                neighbor = &neighbor_block->patches[rpatchno[ir]];
            }
            else
            {
                /* neighbor patch is on a remote processor */
                neighbor = &domain->ghost_patches[rpatchno[ir]];
            }
            neighbor_patches[ir] = neighbor;
        }
        **iface_neighbor_ptr = iface;
        fclaw2d_patch_face_swap(*iface_neighbor_ptr,&rfaceno);
    }
}

/**
 * \ingroup Averaging
 **/
void cb_face_fill(fclaw2d_domain_t *domain,
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

    const amr_options_t *gparms = get_domain_parms(domain);
    const int refratio = gparms->refratio;

    fclaw_bool intersects_phys_bdry[NumFaces];
    fclaw_bool intersects_block[NumFaces];

    fclaw2d_get_physical_bc(domain,this_block_idx,this_patch_idx,
                            intersects_phys_bdry);

    fclaw2d_block_get_block_boundary(domain, this_patch, intersects_block);


    /* Transform data needed at block boundaries */
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;                 /* cell-centered data in this routine. */
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;     /* gets filled in below. */

    ClawPatch *this_cp = fclaw2d_clawpatch_get_cp(this_patch);
    for (int iface = 0; iface < NumFaces; iface++)
    {
        int idir = iface/2;

        fclaw_bool is_block_face;
        fclaw_bool is_interior_face;
        get_face_type(domain,
                      iface,
                      intersects_phys_bdry,
                      intersects_block,
                      &is_block_face,
                      &is_interior_face);


        if (is_interior_face)
        {
            /* Output arguments */
            int neighbor_block_idx;
            int neighbor_level;   /* = -1, 0, 1 */
            int *ref_flag_ptr = &neighbor_level;
            int fine_grid_pos = -1;
            int *fine_grid_pos_ptr = &fine_grid_pos;

            /* Get the face neighbor relative to the neighbor's coordinate
               orientation (this isn't used here) */
            int iface_neighbor;
            int *iface_neighbor_ptr = &iface_neighbor;

            fclaw2d_patch_t* neighbor_patches[p4est_refineFactor];

            /* transform_data.block_iface = iface; */
            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iface, is_block_face,
                               &neighbor_block_idx,
                               neighbor_patches,
                               &ref_flag_ptr,
                               &fine_grid_pos_ptr,
                               &iface_neighbor_ptr,
                               transform_data.transform);

            /* fclaw_bool block_boundary = (neighbor_block_idx >= 0); */
            if (ref_flag_ptr == NULL)
            {
                /* We should never end up here */
                printf("cb_face_fill (fclaw2d_face_neighbors.cpp) : no face found\n");
                exit(0);
            }

            /* Parallel distribution keeps siblings on same processor */
            fclaw_bool remote_neighbor;
            remote_neighbor = fclaw2d_patch_is_ghost(neighbor_patches[0]);
            if (is_coarse)
            {
                if (neighbor_level == FINER_GRID)
                {
                    for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                    {
                        ClawPatch *fine_neighbor_cp =
                            fclaw2d_clawpatch_get_cp(neighbor_patches[igrid]);
                        transform_data.neighbor_patch = neighbor_patches[igrid];
                        transform_data.fine_grid_pos = igrid;
                        if (interpolate_to_neighbor && !remote_neighbor)
                        {
                            /* interpolate to igrid */
                            this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                            refratio,fine_neighbor_cp,
                                                            time_interp,igrid,
                                                            &transform_data);
                        }
                        else if (average_from_neighbor)
                        {
                            /* average from igrid */
                            this_cp->average_face_ghost(idir,iface,p4est_refineFactor,
                                                        refratio,fine_neighbor_cp,
                                                        time_interp,igrid,
                                                        &transform_data);
                        }
                    }
                }
                else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                {
                    /* Copy to same size patch */
                    fclaw2d_patch_t *neighbor_patch = neighbor_patches[0];
                    ClawPatch *neighbor_cp = fclaw2d_clawpatch_get_cp(neighbor_patch);
                    transform_data.neighbor_patch = neighbor_patches[0];
                    transform_data.fine_grid_pos = 0;
                    this_cp->exchange_face_ghost(iface,neighbor_cp,&transform_data);
                }
            }
            else if (is_fine && remote_neighbor && read_parallel_patches)
            {
                /* Swap 'this_patch' and the neighbor patch */
                ClawPatch *coarse_cp = fclaw2d_clawpatch_get_cp(neighbor_patches[0]);
                ClawPatch *fine_cp = this_cp;

                /* Figure out which grid we got */
                int igrid = fine_grid_pos;  /* returned from get_face_neighbors, above */

                int iface_coarse = iface_neighbor;
                int this_face = iface;

                /* make sure we are not using an unintialized variable */
                FCLAW_ASSERT (fine_grid_pos_ptr != NULL);
                FCLAW_ASSERT (fine_grid_pos != -1);

                /* Redo the transformation */
                if (neighbor_block_idx >= 0)
                {
                    fclaw2d_patch_face_transformation (iface_coarse, this_face,
                                                       transform_data.transform);
                }
                transform_data.this_patch = neighbor_patches[0];  /* only one (coarse) neighbor */
                transform_data.neighbor_patch = this_patch;
                transform_data.fine_grid_pos = igrid;

                if (neighbor_level == COARSER_GRID)
                {
                    if (average_from_neighbor)
                    {
                        /* Neighbor patch is a coarser grid;  we want to average 'this_patch' to it */
                        coarse_cp->average_face_ghost(idir,iface_coarse,
                                                      p4est_refineFactor,refratio,
                                                      fine_cp,time_interp,
                                                      igrid, &transform_data);
                    }
                    else if (interpolate_to_neighbor)
                    {
                        /* Neighbor patch is a coarser grid;  we want to interpolate from parallal
                           patch to 'this' patch (the finer grid) */
                        coarse_cp->interpolate_face_ghost(idir,iface_coarse,
                                                          p4est_refineFactor,refratio,
                                                          fine_cp,time_interp,
                                                          igrid, &transform_data);
                    }
                }
                else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                {
                    /* Copy to a parallal patch;  we need these values later for
                       interpolation to on-proc fine grid corners */
                    coarse_cp->exchange_face_ghost(iface_coarse,this_cp,&transform_data);
                }
                else if (neighbor_level == FINER_GRID)
                {
                    /* We don't need to interpolate to parallel patches */
                }
            }
        }  /* End of interior face */
    } /* End of iface loop */
}
