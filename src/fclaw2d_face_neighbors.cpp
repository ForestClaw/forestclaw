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
#include <fclaw2d_patch.h>

#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_clawpatch.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* This is used to determine neighbor patch relative level
   (finer, coarser or samesize) */
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
                        int ftransform[],
                        fclaw2d_transform_data_t* ftransform_finegrid)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    int rproc[p4est_refineFactor];
    int rblockno;
    int rpatchno[p4est_refineFactor];
    int rfaceno;
    int num_neighbors;

    for(int ir = 0; ir < p4est_refineFactor; ir++)
    {
        neighbor_patches[ir] = NULL;
    }

    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_NEIGHBOR_SEARCH]);
    fclaw2d_patch_relation_t neighbor_type =
        fclaw2d_patch_face_neighbors(domain,
                                     this_block_idx,
                                     this_patch_idx,
                                     iface,
                                     rproc,
                                     &rblockno,
                                     rpatchno,
                                     &rfaceno);
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_NEIGHBOR_SEARCH]);


    /* ------------------------------
      neighbor_type is one of :
      FCLAW2D_PATCH_BOUNDARY,    -- physical boundary
      FCLAW2D_PATCH_HALFSIZE,
      FCLAW2D_PATCH_SAMESIZE,
      FCLAW2D_PATCH_DOUBLESIZE
      ------------------------------- */

    if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
    {
        /* This case should be excluded by earlier checks */
        printf("get_face_neighbors (fclaw2d_face_neighbors.cpp) : No patch " \
               "found\n");
        exit(0);
    }
    else
    {
        *neighbor_block_idx = is_block_face ? rblockno : -1;
        /* Get encoding of transforming a neighbor coordinate across a face */
        fclaw2d_patch_face_transformation (iface, rfaceno, ftransform);

        int iface1, rface1;
        iface1 = iface;
        rface1 = rfaceno;
        fclaw2d_patch_face_swap(&iface1,&rface1);
        fclaw2d_patch_face_transformation (iface1, rface1,
                                           ftransform_finegrid->transform);
        ftransform_finegrid->block_iface = iface1;
        **iface_neighbor_ptr = iface1;


        if (!is_block_face)
        {
            /* If we are within one patch this is a special case */
            FCLAW_ASSERT (*neighbor_block_idx == -1);
            fclaw2d_patch_face_transformation_intra (ftransform);
            fclaw2d_patch_face_transformation_intra
                (ftransform_finegrid->transform);
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
    fclaw2d_vtable_t vt = fclaw2d_get_vtable(domain);
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

    fclaw2d_physical_get_bc(domain,this_block_idx,this_patch_idx,
                            intersects_phys_bdry);

    fclaw2d_block_get_block_boundary(domain, this_patch, intersects_block);


    /* Transform data needed at block boundaries */
    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;                 /* cell-centered data in this routine. */
    transform_data.this_patch = this_patch;
    transform_data.neighbor_patch = NULL;     /* gets filled in below. */

    fclaw2d_transform_data_t transform_data_finegrid;
    transform_data_finegrid.mx = gparms->mx;
    transform_data_finegrid.my = gparms->my;
    transform_data_finegrid.based = 1;   /* cell-centered data in this routine. */

#if 0
    ClawPatch *this_cp = fclaw2d_clawpatch_get_cp(this_patch);
#endif
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


        if (is_interior_face)  /* Not on a physical boundary */
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

            /* Reset this in case it got set in a remote copy */
            transform_data.this_patch = this_patch;
#if 0
            this_cp = fclaw2d_clawpatch_get_cp(this_patch);
#endif

            transform_data_finegrid.block_iface = -1;

            /* transform_data.block_iface = iface; */
            get_face_neighbors(domain,
                               this_block_idx,
                               this_patch_idx,
                               iface,
			       is_block_face,
                               &neighbor_block_idx,
                               neighbor_patches,
                               &ref_flag_ptr,
                               &fine_grid_pos_ptr,
                               &iface_neighbor_ptr,
                               transform_data.transform,
                               &transform_data_finegrid);

            /* Needed for switching the context */
            transform_data_finegrid.this_patch = neighbor_patches[0];
            transform_data_finegrid.neighbor_patch = this_patch;


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
                        remote_neighbor = fclaw2d_patch_is_ghost(neighbor_patches[igrid]);
                        if (!((read_parallel_patches && remote_neighbor) || !remote_neighbor))
                        {
                            continue;
                        }
#if 0
                        ClawPatch *fine_neighbor_cp =
                            fclaw2d_clawpatch_get_cp(neighbor_patches[igrid]);
#endif
                        fclaw2d_patch_t* fine_patch = neighbor_patches[igrid];
                        fclaw2d_patch_t *coarse_patch = this_patch;
                        transform_data.neighbor_patch = neighbor_patches[igrid];
                        if (interpolate_to_neighbor && !remote_neighbor)
                        {
                            /* interpolate to igrid */
                            vt.interpolate_face(domain,this_patch,fine_patch,idir,
                                                iface,p4est_refineFactor,refratio,
                                                time_interp,igrid,&transform_data);
#if 0
                            this_cp->interpolate_face_ghost(idir,iface,p4est_refineFactor,
                                                            refratio,fine_neighbor_cp,
                                                            time_interp,igrid,
                                                            &transform_data);
#endif
                        }
                        else if (average_from_neighbor)
                        {
                            /* average from igrid */
                            vt.average_face(domain,coarse_patch,fine_patch,idir,
                                            iface,p4est_refineFactor,
                                            refratio,time_interp,igrid,
                                            &transform_data);
#if 0
                            this_cp->average_face_ghost(idir,iface,p4est_refineFactor,
                                                        refratio,fine_neighbor_cp,
                                                        time_interp,igrid,
                                                        &transform_data);
#endif
                        }
                    }
                }
                else if (neighbor_level == SAMESIZE_GRID && copy_from_neighbor)
                {
                    /* Copy to same size patch */
                    fclaw2d_patch_t *neighbor_patch = neighbor_patches[0];
                    transform_data.neighbor_patch = neighbor_patch;
#if 0
                    ClawPatch *neighbor_cp = fclaw2d_clawpatch_get_cp(neighbor_patch);
                    this_cp->exchange_face_ghost(iface,neighbor_cp,time_interp,
                                                 &transform_data);
#endif
                    vt.copy_face(domain,this_patch,neighbor_patch,iface,
                                        time_interp,&transform_data);

                    /* We also need to copy _to_ the remote neighbor; switch contexts, but
                       use ClawPatches that are only in scope here, to avoid
                       conflicts with above uses of the same variables. This is needed
                       in case we want to interpolate to adjacent corners on fine grids.*/
                    if (remote_neighbor)
                    {
                        /* Create a new transform so we don't mess up the original one */
                        int this_iface = iface_neighbor;
#if 0
                        ClawPatch *neighbor_cp = this_cp;
                        ClawPatch *this_cp = fclaw2d_clawpatch_get_cp(neighbor_patch);
                        this_cp->exchange_face_ghost(this_iface,neighbor_cp,
                                                     time_interp,
                                                     &transform_data_finegrid);
#endif

                        vt.copy_face(domain,neighbor_patch,this_patch,this_iface,
                                            time_interp,&transform_data_finegrid);
                    }
                }
            }
            else if (is_fine && neighbor_level == COARSER_GRID && remote_neighbor
                     && read_parallel_patches)
            {
                int iface_coarse = iface_neighbor;
                int igrid = fine_grid_pos;  /* Not used */
                int idir_coarse = iface_coarse/2;

                /* Swap 'this_patch' (fine grid) and the neighbor patch (a coarse grid) */
#if 0
                ClawPatch *coarse_cp = fclaw2d_clawpatch_get_cp(neighbor_patches[0]);
                ClawPatch *fine_cp = fclaw2d_clawpatch_get_cp(this_patch);
#endif
                fclaw2d_patch_t* coarse_patch = neighbor_patches[0];
                fclaw2d_patch_t* fine_patch = this_patch;

		if (average_from_neighbor)
                {
		    /* Average from 'this' grid (fine grid) to remote grid (coarse grid) */
                    vt.average_face(domain,coarse_patch,fine_patch,idir_coarse,iface_coarse,
                                    p4est_refineFactor,refratio,
                                    time_interp,igrid,&transform_data_finegrid);
#if 0
		    coarse_cp->average_face_ghost(idir_coarse,iface_coarse,
						  p4est_refineFactor,refratio,
						  fine_cp,time_interp,
						  igrid, &transform_data_finegrid);
#endif
                }
                else if (interpolate_to_neighbor)
                {
                    /* Interpolate from remote neighbor to 'this' patch (the finer grid */
                    vt.interpolate_face(domain,coarse_patch,fine_patch,
                                        idir_coarse,iface_coarse,
                                        p4est_refineFactor,refratio,
                                        time_interp,
                                        igrid,&transform_data_finegrid);
#if 0
                    coarse_cp->interpolate_face_ghost(idir_coarse,iface_coarse,
                                                      p4est_refineFactor,refratio,
                                                      fine_cp,time_interp,
                                                      igrid, &transform_data_finegrid);
#endif
                }
            }
        }  /* End of interior face */
    } /* End of iface loop */
}


void fclaw2d_face_neighbor_ghost(fclaw2d_domain_t* domain,
                                 int minlevel,
                                 int maxlevel,
                                 int time_interp)
{
    fclaw2d_vtable_t vt = fclaw2d_get_vtable(domain);
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    const amr_options_t *gparms = get_domain_parms(domain);
    int refratio = gparms->refratio;

    int rproc[2];
    int rpatchno[2];
    int rblockno;
    int rfaceno;

    int min_interp_level = time_interp ? minlevel-1 : minlevel;

    fclaw2d_transform_data_t transform_data;
    transform_data.mx = gparms->mx;
    transform_data.my = gparms->my;
    transform_data.based = 1;      /* cell-centered data in this routine. */

    fclaw2d_domain_indirect_t *ind = ddata->domain_indirect;

    /* Loop over ghost patches to do any face exchanges that didn't happen
       before the ghost patch exchange was done */

    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* this_ghost_patch = &domain->ghost_patches[i];
        int level = this_ghost_patch->level;
        if (level < min_interp_level)
        {
            /* We don't need to worry about ghost patches that are at
               coarser levels than we are currently working on */
            continue;
        }

        int use_timeinterp_patch = 0;
        if (time_interp && level == minlevel-1)
        {
            use_timeinterp_patch = 1;
        }


#if 0
        ClawPatch *this_cp = fclaw2d_clawpatch_get_cp(this_ghost_patch);
#endif

        int this_ghost_idx = i;

        transform_data.this_patch = this_ghost_patch;

        for (int iface = 0; iface < NumFaces; iface++)
        {
            int idir = iface/2;

            /* We are only looking for faces between two ghost patches
               from different processors, since these will not have
               exchange face data before being thrown over proc fence.
            */
            fclaw2d_patch_relation_t neighbor_type =
                fclaw2d_domain_indirect_neighbors(domain,
                                                  ind,
                                                  this_ghost_idx,
                                                  iface,rproc,
                                                  &rblockno, rpatchno,
                                                  &rfaceno);

            if (neighbor_type != FCLAW2D_PATCH_BOUNDARY)
            {
                /* We have a neighbor ghost patch that came from a
                   different proc */

                fclaw_bool intersects_block[NumFaces];
                fclaw2d_block_get_block_boundary(domain, this_ghost_patch,
                                                 intersects_block);
                int is_block_face = intersects_block[iface];


                fclaw2d_patch_face_transformation (iface, rfaceno,
                                                   transform_data.transform);

                if (!is_block_face)
                {
                    fclaw2d_patch_face_transformation_intra(transform_data.transform);
                }
                if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
                {
                    /* Copy from same size neighbor */
                    fclaw2d_patch_t *neighbor_patch = &domain->ghost_patches[rpatchno[0]];
                    transform_data.neighbor_patch = neighbor_patch;
#if 0
                    ClawPatch *neighbor_cp = fclaw2d_clawpatch_get_cp(neighbor_patch);
                    this_cp->exchange_face_ghost(iface,neighbor_cp,
                                                 use_timeinterp_patch,
                                                 &transform_data);
#endif
                    vt.copy_face(domain,this_ghost_patch,neighbor_patch,iface,
                                        use_timeinterp_patch,&transform_data);

                    ++ddata->count_multiproc_corner;
                }
                else if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
                {
                    /* Average from fine grid neighbor */
                    for (int igrid = 0; igrid < p4est_refineFactor; igrid++)
                    {
                        if (rpatchno[igrid] != -1)
                        {
                            fclaw2d_patch_t* coarse_patch = this_ghost_patch;
                            fclaw2d_patch_t *fine_patch =
                                &domain->ghost_patches[rpatchno[igrid]];
                            transform_data.neighbor_patch = fine_patch;
                            vt.average_face(domain,coarse_patch,fine_patch,
                                            idir,iface,p4est_refineFactor,
                                            refratio,use_timeinterp_patch,
                                            igrid,&transform_data);
#if 0
                            ClawPatch *neighbor_cp =
                                fclaw2d_clawpatch_get_cp(neighbor_patch);

                            /* Note : igrid isn't used here; could be removed */
                            this_cp->average_face_ghost(idir,iface,p4est_refineFactor,
                                                        refratio,neighbor_cp,
                                                        use_timeinterp_patch,
                                                        igrid,
                                                        &transform_data);
#endif
                        }
                    }
                    ++ddata->count_multiproc_corner;
                }
                else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
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

        }
    }
}

static
void cb_set_neighbor_types(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int blockno,
                           int patchno,
                           void *user)
{

    for (int iface = 0; iface < 4; iface++)
    {
        int rproc[2];
        int rblockno;
        int rpatchno[2];
        int rfaceno;

        fclaw2d_patch_relation_t neighbor_type =
            fclaw2d_patch_face_neighbors(domain,
                                         blockno,
                                         patchno,
                                         iface,
                                         rproc,
                                         &rblockno,
                                         rpatchno,
                                         &rfaceno);

        fclaw2d_patch_set_face_type(this_patch,iface,neighbor_type);
    }

    for (int icorner = 0; icorner < 4; icorner++)
    {
        int rproc_corner;
        int cornerpatchno;
        int cornerblockno;
        int rcornerno;
        fclaw2d_patch_relation_t neighbor_type;

        int has_corner_neighbor =
            fclaw2d_patch_corner_neighbors(domain,
                                           blockno,
                                           patchno,
                                           icorner,
                                           &rproc_corner,
                                           &cornerblockno,
                                           &cornerpatchno,
                                           &rcornerno,
                                           &neighbor_type);

        fclaw2d_patch_set_corner_type(this_patch,icorner,neighbor_type);
        if (!has_corner_neighbor)
        {
            fclaw2d_patch_set_missing_corner(this_patch,icorner);
        }
    }
    fclaw2d_patch_neighbors_set(this_patch);
}

void fclaw2d_regrid_set_neighbor_types(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_NEIGHBOR_SEARCH]);
    fclaw2d_domain_iterate_patches(domain,cb_set_neighbor_types,(void*) NULL);
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_NEIGHBOR_SEARCH]);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
