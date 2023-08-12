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

/** \file fclaw2d_face_neighbors.cpp
 * Average, coarsen and copy between grids at faces.
 **/

#include <fclaw2d_face_neighbors.h>

#include <fclaw_global.h>
#include <fclaw2d_defs.h>

#include <fclaw2d_block.h>
#include <fclaw_patch.h>
#include <fclaw_ghost_fill.h>
#include <fclaw_options.h>
#include <fclaw_physical_bc.h>
#include <fclaw_regrid.h>
#include <fclaw_domain.h>


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
void get_face_type(fclaw_global_t *glob,
				   int iface,
				   int intersects_phys_bdry[],
				   int intersects_block[],
				   int *is_block_face,
				   int *is_interior_face)
{
	*is_block_face = intersects_block[iface];
	*is_interior_face = !intersects_phys_bdry[iface];
}


static
void get_face_neighbors(fclaw_global_t *glob,
						int this_block_idx,
						int this_patch_idx,
						int iface,
						int is_block_face,
						int *neighbor_block_idx,
						fclaw_patch_t* neighbor_patches[],						
						int **ref_flag_ptr,
						int **fine_grid_pos_ptr,
						int **iface_neighbor_ptr,
						int ftransform[],
						fclaw_patch_transform_data_t* ftransform_finegrid)
{
	fclaw_domain_t *domain = glob->domain;
	int rproc[FCLAW2D_NUMFACENEIGHBORS];
	int rblockno;
	int rpatchno[FCLAW2D_NUMFACENEIGHBORS];
	int rfaceno;
	int num_neighbors;
	int ir;

	for(ir = 0; ir < FCLAW2D_NUMFACENEIGHBORS; ir++)
	{
		neighbor_patches[ir] = NULL;
	}

	fclaw_timer_start (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);
	fclaw_patch_relation_t neighbor_type =
	fclaw_patch_face_neighbors(domain,
								 this_block_idx,
								 this_patch_idx,
								 iface,
								 rproc,
								 &rblockno,
								 rpatchno,
								 &rfaceno);
	fclaw_timer_stop (&glob->timers[FCLAW_TIMER_NEIGHBOR_SEARCH]);


	/* ------------------------------
	  neighbor_type is one of :
	  FCLAW2D_PATCH_BOUNDARY,    -- physical boundary
	  FCLAW2D_PATCH_HALFSIZE,
	  FCLAW2D_PATCH_SAMESIZE,
	  FCLAW2D_PATCH_DOUBLESIZE
	  ------------------------------- */

	if (neighbor_type == FCLAW_PATCH_BOUNDARY)
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
		fclaw_patch_transform_blockface (glob, iface, rfaceno, ftransform);

		int iface1, rface1;
		iface1 = iface;
		rface1 = rfaceno;
		fclaw2d_patch_face_swap(&iface1,&rface1);
		fclaw_patch_transform_blockface (glob, iface1, rface1,
										   ftransform_finegrid->transform);
		ftransform_finegrid->d2->block_iface = iface1;
		**iface_neighbor_ptr = iface1;


		if (!is_block_face)
		{
			/* If we are within one patch this is a special case */
			FCLAW_ASSERT (*neighbor_block_idx == -1);
			fclaw_patch_transform_blockface_intra (glob, ftransform);
			fclaw_patch_transform_blockface_intra
				(glob, ftransform_finegrid->transform);
		}

		if (neighbor_type == FCLAW_PATCH_SAMESIZE)
		{
			**ref_flag_ptr = 0;
			*fine_grid_pos_ptr = NULL;
			num_neighbors = 1;
		}
		else if (neighbor_type == FCLAW_PATCH_DOUBLESIZE)
		{
			**ref_flag_ptr = -1;
			**fine_grid_pos_ptr = rproc[1];    /* Special storage for fine grid info */
			num_neighbors = 1;
		}
		else if (neighbor_type == FCLAW_PATCH_HALFSIZE)
		{
			/* Patch has two neighbors */
			**ref_flag_ptr = 1; /* patches are at one level finer */
			*fine_grid_pos_ptr = NULL;
			num_neighbors = FCLAW2D_NUMFACENEIGHBORS;
		}
		else
		{
			printf ("Illegal fclaw2d_patch_face_neighbors return value\n");
			exit (1);
		}

		for(ir = 0; ir < num_neighbors; ir++)
		{
			fclaw_patch_t *neighbor;
			if (rproc[ir] == domain->mpirank)
			{
				/* neighbor patch is local */
				fclaw_block_t *neighbor_block = &domain->blocks[rblockno];
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
void cb_face_fill(fclaw_domain_t *domain,
				  fclaw_patch_t *this_patch,
				  int this_block_idx,
				  int this_patch_idx,
				  void *user)
{
	fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user; 

	fclaw_exchange_info_t *filltype = (fclaw_exchange_info_t*) s->user;
	int time_interp = filltype->time_interp;
	int is_coarse = filltype->grid_type == FCLAW_IS_COARSE;
	int is_fine = filltype->grid_type == FCLAW_IS_FINE;

	int iface, igrid;

	int read_parallel_patches = filltype->read_parallel_patches;

	int copy_from_neighbor = filltype->exchange_type == FCLAW_COPY;
	int average_from_neighbor = filltype->exchange_type == FCLAW_AVERAGE;
	int interpolate_to_neighbor = filltype->exchange_type == FCLAW_INTERPOLATE;

	int time_sync_fine_to_coarse = filltype->exchange_type == FCLAW_TIME_SYNC_F2C;
	int time_sync_samesize = filltype->exchange_type == FCLAW_TIME_SYNC_SAMESIZE;

	const fclaw_options_t *gparms = fclaw_get_options(s->glob);
	const int refratio = gparms->refratio;

	int intersects_phys_bdry[FCLAW2D_NUMFACES];
	int intersects_block[FCLAW2D_NUMFACES];

	fclaw_physical_get_bc(s->glob,this_block_idx,this_patch_idx,
							intersects_phys_bdry);

	fclaw2d_block_get_block_boundary(s->glob, this_patch, intersects_block);


	/* Transform data needed at block boundaries */
	fclaw_patch_transform_data_t transform_data;

	transform_data.dim = s->glob->domain->dim;

	fclaw_patch_transform_data_d2_t transform_data_d2;
	transform_data.d2 = &transform_data_d2;
#ifndef P4_TO_P8
	transform_data.d3 = NULL;
#else
	transform_data.d2 = NULL;
#endif

	transform_data.glob = s->glob;
	transform_data.based = 1;             /* Set by user defined patch routine */
	transform_data.this_patch = this_patch;
	transform_data.neighbor_patch = NULL;     /* gets filled in below. */

	/* This calls a patch-specific initialization routine  - does nothing yet. */
	fclaw_patch_transform_init_data(s->glob,this_patch,
									  this_block_idx,
									  this_patch_idx,
									  &transform_data);

	fclaw_patch_transform_data_t transform_data_finegrid;

	transform_data_finegrid.dim = s->glob->domain->dim;

	fclaw_patch_transform_data_d2_t transform_data_finegrid_d2;
	transform_data_finegrid.d2 = &transform_data_finegrid_d2;
#ifndef P4_TO_P8
	transform_data_finegrid.d3 = NULL;
#else
	transform_data_finegrid.d2 = NULL;
#endif

	transform_data_finegrid.glob = s->glob;
	transform_data_finegrid.based = 1;   /* cell-centered data in this routine. */

	/* This calls a patch-specific initialization routine - does nothing yet. */
	fclaw_patch_transform_init_data(s->glob,this_patch,
									  this_block_idx,
									  this_patch_idx,
									  &transform_data_finegrid);

	for (iface = 0; iface < FCLAW2D_NUMFACES; iface++)
	{
		int idir = iface/2;

		int is_block_face;
		int is_interior_face;
		get_face_type(s->glob,
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

			fclaw_patch_t* neighbor_patches[FCLAW2D_NUMFACENEIGHBORS];

			/* Reset this in case it got set in a remote copy */
			transform_data.this_patch = this_patch;
			transform_data_finegrid.d2->block_iface = -1;

			/* transform_data.block_iface = iface; */
			get_face_neighbors(s->glob,
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


			/* int block_boundary = (neighbor_block_idx >= 0); */
			if (ref_flag_ptr == NULL)
			{
				/* We should never end up here */
				printf("cb_face_fill (fclaw2d_face_neighbors.cpp) : no face found\n");
				exit(0);
			}

			/* Parallel distribution keeps siblings on same processor */
			int remote_neighbor;
			remote_neighbor = fclaw_patch_is_ghost(neighbor_patches[0]);
			if (is_coarse)
			{
				if (neighbor_level == FINER_GRID)
				{
					for (igrid = 0; igrid < 2; igrid++)
					{
						remote_neighbor = fclaw_patch_is_ghost(neighbor_patches[igrid]);
						int valid_remote = read_parallel_patches && remote_neighbor;
						int local_neighbor = !remote_neighbor;
						if (!(local_neighbor || valid_remote))
						{
							continue;
						}
						fclaw_patch_t* fine_patch = neighbor_patches[igrid];
						fclaw_patch_t *coarse_patch = this_patch;
						transform_data.neighbor_patch = neighbor_patches[igrid];
						if (time_sync_fine_to_coarse)
						{
							int coarse_patchno = this_patch_idx;
							int coarse_blockno = this_block_idx;
							int fine_blockno = neighbor_block_idx;
							
							/* Add correction to coarse grid from fine grid.  */
							fclaw_patch_time_sync_f2c(s->glob,coarse_patch,fine_patch,
							                            coarse_blockno, fine_blockno, 
							                            coarse_patchno,
							                            idir,igrid,
							                            iface,time_interp,
							                            &transform_data);
						}
						else if (interpolate_to_neighbor && !remote_neighbor)
						{
							/* interpolate to igrid */
							fclaw_patch_interpolate_face(s->glob,this_patch,fine_patch,
														   idir,iface,FCLAW2D_REFINEFACTOR,
														   refratio,time_interp,igrid,
														   &transform_data);
						}
						else if (average_from_neighbor)
						{
							/* average from igrid */
							fclaw_patch_average_face(s->glob,coarse_patch,fine_patch,idir,
													   iface,FCLAW2D_REFINEFACTOR,
													   refratio,time_interp,igrid,
													   &transform_data);
						}
					}
				}
				else if (neighbor_level == SAMESIZE_GRID && 
				         (copy_from_neighbor || time_sync_samesize))
				{
					/* Copy to same size patch */
					fclaw_patch_t *neighbor_patch = neighbor_patches[0];
					transform_data.neighbor_patch = neighbor_patch;

					if (time_sync_samesize)    // && is_block_face)
					{
						/* Correct for metric discontinuities at block boundaries */
						if (!fclaw_patch_is_ghost(this_patch))
						{							
							fclaw_patch_time_sync_samesize(s->glob,this_patch,
							                                 neighbor_patch,
							                                 iface,idir,
							                                 &transform_data);
						}
					}
					else
					{                        
						fclaw_patch_copy_face(s->glob,this_patch,neighbor_patch,iface,
												time_interp,&transform_data);
					}

					/* We also need to copy _to_ the remote neighbor; switch contexts, but
					   use ClawPatches that are only in scope here, to avoid
					   conflicts with above uses of the same variables. This is needed
					   in case we want to interpolate to adjacent corners on fine grids.*/
					if (remote_neighbor)
					{
						/* Create a new transform so we don't mess up the original one */
						int this_iface = iface_neighbor;

						if (time_sync_samesize) // && is_block_face)
						{
							/* Correct ghost patches, since these will be used to copy or
							interpolate to local grids. */
							fclaw_patch_time_sync_samesize(s->glob,neighbor_patch,this_patch,
							                                 this_iface,idir,
							                                 &transform_data_finegrid);                            
						}
						else
						{
							fclaw_patch_copy_face(s->glob,neighbor_patch,this_patch,
													this_iface, time_interp,
													&transform_data_finegrid);                            
						}
					}
				}
			}
			else if (is_fine && neighbor_level == COARSER_GRID && remote_neighbor
					 && read_parallel_patches)
			{
				int iface_coarse = iface_neighbor;
				int igrid = fine_grid_pos;  /* Not used */
				int idir_coarse = iface_coarse/2;

				/* Swap 'this_patch' (fine grid) and the neighbor patch 
				(a coarse grid) */
				fclaw_patch_t* coarse_patch = neighbor_patches[0];
				fclaw_patch_t* fine_patch = this_patch;

				if (average_from_neighbor)
				{
						/* Average from 'this' grid (fine grid) to remote grid 
					(coarse grid) */
					fclaw_patch_average_face(s->glob,coarse_patch,fine_patch,
											   idir_coarse,iface_coarse,
											   FCLAW2D_REFINEFACTOR,refratio,
											   time_interp,igrid,
											   &transform_data_finegrid);
				}
				else if (interpolate_to_neighbor)
				{
					/* Interpolate from remote neighbor to 'this' patch (the finer grid */
					fclaw_patch_interpolate_face(s->glob,coarse_patch,fine_patch,
												   idir_coarse,iface_coarse,
												   FCLAW2D_REFINEFACTOR,refratio,
												   time_interp,
												   igrid,&transform_data_finegrid);
				}
			}
		}  /* End of interior face */
	} /* End of iface loop */
}


void fclaw_face_neighbor_ghost(fclaw_global_t* glob,
								 int minlevel,
								 int maxlevel,
								 int time_interp)
{
	fclaw_domain_t *domain = glob->domain;

	fclaw2d_domain_data_t *ddata = domain->d2;
	const fclaw_options_t *gparms = fclaw_get_options(glob);
	int refratio = gparms->refratio;

	int rproc[2];
	int rpatchno[2];
	int rblockno;
	int rfaceno;
	int i, iface, igrid;

	int min_interp_level = time_interp ? minlevel-1 : minlevel;

	fclaw_patch_transform_data_t transform_data;

	transform_data.dim = domain->dim;

	fclaw_patch_transform_data_d2_t transform_data_d2;
	transform_data.d2 = &transform_data_d2;
#ifndef P4_TO_P8
	transform_data.d3 = NULL;
#else
	transform_data.d2 = NULL;
#endif

	transform_data.glob = glob;
	transform_data.based = 1;      /* cell-centered data in this routine. */

	fclaw2d_domain_indirect_t *ind = ddata->domain_indirect;

	/* Loop over ghost patches to do any face exchanges that didn't happen
	   before the ghost patch exchange was done */

	for(i = 0; i < domain->num_ghost_patches; i++)
	{
		fclaw_patch_t* this_ghost_patch = &domain->ghost_patches[i];
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


		int this_ghost_idx = i;

		transform_data.this_patch = this_ghost_patch;
		fclaw_patch_transform_init_data(glob,this_ghost_patch,
		                              	  -1,
		                                  i,
		                                  &transform_data);

		for (iface = 0; iface < FCLAW2D_NUMFACES; iface++)
		{
			int idir = iface/2;

			/* We are only looking for faces between two ghost patches
			   from different processors, since these will not have
			   exchange face data before being thrown over proc fence.
			*/
			fclaw_patch_relation_t neighbor_type =
				fclaw2d_domain_indirect_neighbors(domain,
												  ind,
												  this_ghost_idx,
												  iface,rproc,
												  &rblockno, rpatchno,
												  &rfaceno);

			if (neighbor_type != FCLAW_PATCH_BOUNDARY)
			{
				/* We have a neighbor ghost patch that came from a
				   different proc */

				int intersects_block[FCLAW2D_NUMFACES];
				fclaw2d_block_get_block_boundary(glob, this_ghost_patch,
												 intersects_block);
				int is_block_face = intersects_block[iface];


				fclaw_patch_transform_blockface (glob, iface, rfaceno,
												   transform_data.transform);

				if (!is_block_face)
				{
					fclaw_patch_transform_blockface_intra(glob, transform_data.transform);
				}
				if (neighbor_type == FCLAW_PATCH_SAMESIZE)
				{
					/* Copy from same size neighbor */
					fclaw_patch_t *neighbor_patch = &domain->ghost_patches[rpatchno[0]];
					transform_data.neighbor_patch = neighbor_patch;
					fclaw_patch_copy_face(glob,this_ghost_patch,neighbor_patch,iface,
											use_timeinterp_patch,&transform_data);

					++glob->count_multiproc_corner;
				}
				else if (neighbor_type == FCLAW_PATCH_HALFSIZE)
				{
					/* Average from fine grid neighbor */
					for (igrid = 0; igrid < FCLAW2D_NUMFACENEIGHBORS; igrid++)
					{
						if (rpatchno[igrid] != -1)
						{
							fclaw_patch_t* coarse_patch = this_ghost_patch;
							fclaw_patch_t *fine_patch =
								&domain->ghost_patches[rpatchno[igrid]];
							transform_data.neighbor_patch = fine_patch;
							fclaw_patch_average_face(glob,coarse_patch,fine_patch,
													   idir,iface,FCLAW2D_REFINEFACTOR,
													   refratio,use_timeinterp_patch,
													   igrid,&transform_data);
						}
					}
					++glob->count_multiproc_corner;
				}
				else if (neighbor_type == FCLAW_PATCH_DOUBLESIZE)
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
