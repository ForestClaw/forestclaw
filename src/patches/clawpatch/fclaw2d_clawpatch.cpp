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

#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 2
#endif

#if REFINE_DIM == 2 && PATCH_DIM == 2

#include <fclaw2d_clawpatch.h>
#include <fclaw_clawpatch.hpp>

#include <fclaw_clawpatch_diagnostics.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw2d_clawpatch_fort.h>
#include <fclaw2d_clawpatch_conservation.h>
#include <fclaw2d_clawpatch_conservation_fort.h>
#include <fclaw2d_clawpatch_transform.h>
#include <fclaw_clawpatch_pillow.h>  

#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>

#include <fclaw2d_metric.h>
#include <fclaw2d_metric.hpp>

#elif REFINE_DIM == 2 && PATCH_DIM == 3

#include <fclaw3dx_clawpatch.h>
#include <fclaw_clawpatch.hpp>

#include <fclaw_clawpatch_diagnostics.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw3dx_clawpatch_fort.h>
#include <fclaw3dx_clawpatch_transform.h>
#include <fclaw_clawpatch_pillow.h>  

#include <fclaw3dx_clawpatch46_fort.h>

#include <fclaw3d_metric.h>
#include <fclaw3d_metric.hpp>

#include <_fclaw2d_to_fclaw3dx.h>
#include <_fclaw2d_to_fclaw3d.h>

#endif

#include <fclaw2d_patch.h>  /* Needed to get enum for build modes */

#include <fclaw2d_defs.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_options.h>

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_diagnostics.h>

#include <fclaw2d_map_query.h>

#include <fclaw_pointer_map.h>



/* ------------------------------- Static function defs ------------------------------- */

/* Added to turn off time_interp */
static int fill_ghost(fclaw2d_global_t* glob, int time_interp)
{
	const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
	if (fclaw_opt->timeinterp2fillghost)
		/* This will always fill ghost cells using data from "qsync", which is either 
		   coarse grid data at time step, or time interpolated data */
		return 1;
	else
		/* Only fill ghost cells with neighboring data if not doing time interpolation. 
		  If doing time interpolation, then fine grid data at intermediate time steps 
		  will be filled in some other way (.e.g. by advancing the solution in ghost 
		  cells.) */
		return !time_interp;
}


/* Store virtual table for retrieval from anywhere */

static
fclaw_clawpatch_t* get_clawpatch(fclaw2d_patch_t *patch)
{
	fclaw_clawpatch_t *cp = (fclaw_clawpatch_t*) 
					 fclaw2d_patch_get_user_patch(patch);
	return cp;
}

/* Needed for virtual patch function so that the metric class can be independent 
of a clawpatch object */
static 
void* clawpatch_get_metric_patch(fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
#   if PATCH_DIM == 2
		return cp->d2->mp;
#   else
		return cp->d3->mp;
#   endif
}


static
fclaw2d_metric_patch_t* get_metric_patch(fclaw2d_patch_t *patch)
{
	return (fclaw2d_metric_patch_t*) clawpatch_get_metric_patch(patch);
}

/* Return a pointer to either time interpolated data or regular grid data */
static 
double* q_time_sync(fclaw2d_patch_t* patch, int time_interp)
{
	fclaw_clawpatch_t* cp = get_clawpatch(patch);
	if (time_interp)
		return cp->griddata_time_interpolated.dataPtr();
	else
		return cp->griddata.dataPtr();
}

static 
double* clawpatch_get_area(struct fclaw2d_global* glob,
	                       fclaw2d_patch_t* patch)
{
	return fclaw2d_metric_patch_get_area(glob, patch);
}


/* ----------------------------- Creating/deleting patches ---------------------------- */

static
void* clawpatch_new()
{
	fclaw_clawpatch_t *cp = new fclaw_clawpatch_t;    
#if PATCH_DIM == 2
	cp->dim = 2;
	cp->d2 = new fclaw_clawpatch_2d_t;
	/* This patch will only be defined if we are on a manifold. */
	cp->d2->mp = fclaw2d_metric_patch_new();
#else
	cp->dim = 3;
	cp->d3 = new fclaw_clawpatch_3d_t;
	/* This patch will only be defined if we are on a manifold. */
	cp->d3->mp = fclaw2d_metric_patch_new();
#endif

	return (void*) cp;
}

static
void clawpatch_delete(void *patchcp)
{
	FCLAW_ASSERT(patchcp != NULL);
	fclaw_clawpatch_t* cp = (fclaw_clawpatch_t*) patchcp;

#if PATCH_DIM == 2
	fclaw2d_clawpatch_time_sync_delete(&cp->d2->registers);

	FCLAW_ASSERT(cp->d2->mp != NULL);
	fclaw2d_metric_patch_delete(&cp->d2->mp);

	delete cp->d2;
#else
	FCLAW_ASSERT(cp->d3->mp != NULL);
	fclaw2d_metric_patch_delete(&cp->d3->mp);

	delete cp->d3;
#endif

	delete cp;
	patchcp = NULL;
}

/* Maybe this should just be a 'build' function? */
static
void clawpatch_define(fclaw2d_global_t* glob,
					  fclaw2d_patch_t *patch,
					  int blockno, int patchno,
					  fclaw2d_build_mode_t build_mode)
{
	/* We are getting closer to getting rid the class fclaw_clawpatch_t */
	fclaw_clawpatch_t *cp = get_clawpatch(patch);

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	const fclaw_clawpatch_options_t *clawpatch_opt = 
	                     fclaw_clawpatch_get_options(glob);

#   if PATCH_DIM == 2
    	cp->d2->mx = clawpatch_opt->d2->mx;
    	cp->d2->my = clawpatch_opt->d2->my;
#   else
    	cp->d3->mx = clawpatch_opt->d3->mx;
    	cp->d3->my = clawpatch_opt->d3->my;
		cp->d3->mz = clawpatch_opt->d3->mz;
#   endif

	cp->mbc = clawpatch_opt->mbc;
	cp->blockno = blockno;
	cp->meqn = clawpatch_opt->meqn;
	cp->maux = clawpatch_opt->maux;
	cp->mfields = clawpatch_opt->rhs_fields;

	for (int icorner=0; icorner < 4; icorner++)
	{
		fclaw2d_patch_set_block_corner_count(glob,patch,icorner,0);
	}

	fclaw2d_map_context_t* cont = glob->cont;

	int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

	cp->manifold = fclaw_opt->manifold;
	if (cp->manifold)
	{
#       if PATCH_DIM == 2
			cp->d2->xlower = patch->xlower;
			cp->d2->ylower = patch->ylower;
			cp->d2->xupper = patch->xupper;
			cp->d2->yupper = patch->yupper;
#       else
			cp->d3->xlower = patch->xlower;
			cp->d3->ylower = patch->ylower;
			cp->d3->xupper = patch->xupper;
			cp->d3->yupper = patch->yupper;
#           if REFINE_DIM == 2
		    	cp->d3->zlower = 0;
		    	cp->d3->zupper = 1;
#			else
				fclaw_global_essentialf("clawpatch::define : Octrees not yet implemented.\n");
				exit(1);
#      		endif
#      endif

	}
	else
	{
		double ax = fclaw_opt->ax;
		double bx = fclaw_opt->bx;
		double ay = fclaw_opt->ay;
		double by = fclaw_opt->by;

		double xl = patch->xlower;
		double yl = patch->ylower;
		double xu = patch->xupper;
		double yu = patch->yupper;

		double xlower, ylower, xupper, yupper;

		if (is_brick)
		{
			double z;
			/* Scale to [0,1]x[0,1], based on blockno */
			fclaw2d_map_c2m_nomap_brick(cont,cp->blockno,xl,yl,&xlower,&ylower,&z);
			fclaw2d_map_c2m_nomap_brick(cont,cp->blockno,xu,yu,&xupper,&yupper,&z);
		}
		else
		{
			xlower = xl;
			ylower = yl;
			xupper = xu;
			yupper = yu;
		}
#		if PATCH_DIM == 2
			cp->d2->xlower = ax + (bx - ax)*xlower;
			cp->d2->xupper = ax + (bx - ax)*xupper;
			cp->d2->ylower = ay + (by - ay)*ylower;
			cp->d2->yupper = ay + (by - ay)*yupper;
#		else
			cp->d3->xlower = ax + (bx - ax)*xlower;
			cp->d3->xupper = ax + (bx - ax)*xupper;
			cp->d3->ylower = ay + (by - ay)*ylower;
			cp->d3->yupper = ay + (by - ay)*yupper;
			/* Use [az,bz] to scale in z direction.  This should work for 
		    both extruded mesh and octree mesh. */
		    double az = fclaw_opt->az;
		    double bz = fclaw_opt->bz;
			double zlower = 0;		
			double zupper = 1;
			cp->d3->zlower = az + (bz - az)*zlower;
			cp->d3->zupper = az + (bz - az)*zupper;
#		endif

#if REFINE_DIM == 3
#error "clawpatch::define : Octrees not yet implemented."
#endif		

	}

#	if PATCH_DIM == 2
		cp->d2->dx = (cp->d2->xupper - cp->d2->xlower)/cp->d2->mx;
		cp->d2->dy = (cp->d2->yupper - cp->d2->ylower)/cp->d2->my;
# 	else
		cp->d3->dx = (cp->d3->xupper - cp->d3->xlower)/cp->d3->mx;
		cp->d3->dy = (cp->d3->yupper - cp->d3->ylower)/cp->d3->my;
		cp->d3->dz = (cp->d3->zupper - cp->d3->zlower)/cp->d3->mz;
#	endif



	int ll[PATCH_DIM];
	int ur[PATCH_DIM];
	for (int idir = 0; idir < PATCH_DIM; idir++)
	{
		ll[idir] = 1-cp->mbc;
	}
#	if PATCH_DIM == 2
		ur[0] = cp->d2->mx + cp->mbc;
		ur[1] = cp->d2->my + cp->mbc;
#	else
		ur[0] = cp->d3->mx + cp->mbc;
		ur[1] = cp->d3->my + cp->mbc;
		ur[2] = cp->d3->mz + cp->mbc;
#   endif

	Box box(ll,ur,PATCH_DIM);	

	// This will destroy any existing memory n griddata.
	cp->griddata.define(box, cp->meqn);
	if (fclaw_opt->subcycle)
		cp->griddata_time_interpolated.define(box, cp->meqn);

	if (fclaw_opt->compute_error)
	{
		cp->griderror.define(box,cp->meqn);
		cp->exactsolution.define(box,cp->meqn);
	}

	if (clawpatch_opt->maux > 0)
	{		
		cp->aux.define(box,cp->maux);
		if (clawpatch_opt->save_aux)
			cp->aux_save.define(box,cp->maux);
	}

	if (clawpatch_opt->rhs_fields > 0)
	{
		cp->rhs.define(box,cp->mfields);
		if (fclaw_opt->compute_error)
		{
			cp->elliptic_error.define(box,cp->mfields);
			cp->elliptic_soln.define(box,cp->mfields);
		}
	}

	if (fclaw_opt->manifold)
	{
		/*  
			Allocate space for all metric terms. 

		    Note: We pass in detailed info so that the metric patch 
		    doesn't have to know about a clawpatch 
		*/

#if PATCH_DIM == 2
		fclaw2d_metric_patch_define(glob,patch, 
		                            cp->d2->mx, cp->d2->my, cp->mbc, 
		                            cp->d2->dx, cp->d2->dy,
		                            cp->d2->xlower,cp->d2->ylower, cp->d2->xupper, cp->d2->yupper, 
		                            blockno, patchno, build_mode);
#elif PATCH_DIM == 3
		fclaw3d_metric_patch_define(glob,patch, 
		                            cp->d3->mx, cp->d3->my, cp->d3->mz, cp->mbc, 
		                            cp->d3->dx, cp->d3->dy, cp->d3->dz,
		                            cp->d3->xlower,cp->d3->ylower, cp->d3->zlower, 
		                            cp->d3->xupper, cp->d3->yupper, cp->d3->zupper,
		                            blockno, patchno, build_mode);
#endif
	}
	
	/* Build interface registers needed for conservation */
#	if PATCH_DIM == 2
		fclaw2d_clawpatch_time_sync_new(glob,patch,
									  blockno,patchno,&cp->d2->registers);
#	endif

	if (build_mode != FCLAW2D_BUILD_FOR_UPDATE)
		/* If we are building ghost patches, we don't need all the patch memory */
		return;

	cp->griddata_last.define(box, cp->meqn);
	cp->griddata_save.define(box, cp->meqn);

}

static
void clawpatch_build(fclaw2d_global_t *glob,
					 fclaw2d_patch_t *patch,
					 int blockno,
					 int patchno,
					 void *user)
{
	fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

	/* Allocate space for all required arrays */
	clawpatch_define(glob,patch,blockno,patchno,build_mode);

	if (fclaw_opt->manifold)
	{ 
		/* Computes all metric terms : (xp,yp,zp,xd,yd,zd), 
		areas/volumes, basis vectors and (in 2d) curvatures
		and surface normals. 

		Note : Area/volume computations are expensive, since 
		they are all computed on finest levels 

		Note : When building from fine, call routine below 
		where areas/volumes are averaged from fine patches */
		fclaw2d_metric_patch_build(glob,patch,blockno,patchno);
	}

#if PATCH_DIM == 2
	/* Setup for conservation correction */
	fclaw2d_clawpatch_time_sync_setup(glob,patch,blockno,patchno);
#endif
}

static
void clawpatch_build_from_fine(fclaw2d_global_t *glob,
							   fclaw2d_patch_t *fine_patches,
							   fclaw2d_patch_t *coarse_patch,
							   int blockno,
							   int coarse_patchno,
							   int fine0_patchno,
							   fclaw2d_build_mode_t build_mode)
{
	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

	clawpatch_define(glob,coarse_patch,blockno,coarse_patchno,build_mode);

	if (fclaw_opt->manifold)
	{
		/* Average areas/volume from fine grid to coarse grid */
		fclaw2d_metric_patch_build_from_fine(glob, fine_patches, 
		                                     coarse_patch,
											 blockno, coarse_patchno, 
											 fine0_patchno);
	}

#if PATCH_DIM == 2
	/* Setup for conservation correction */
	fclaw2d_clawpatch_time_sync_setup(glob,coarse_patch,blockno,coarse_patchno);
#endif
}


/* -------------------------------- time stepping ------------------------------------- */

static
void clawpatch_save_step(fclaw2d_global_t* glob,
						 fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	cp->griddata_save = cp->griddata;

	/* Some aux arrays are time dependent, or contain part of the solution.  In this case, 
	   we should save the aux array in case we need to re-take a time step */
	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
	if (clawpatch_opt->save_aux)
		cp->aux_save = cp->aux;
}


static
void clawpatch_restore_step(fclaw2d_global_t* glob,
							fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	cp->griddata = cp->griddata_save;

	/* Restore the aux array after before retaking a time step */
	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
	if (clawpatch_opt->save_aux)
		cp->aux = cp->aux_save;
}

static
void clawpatch_setup_timeinterp(fclaw2d_global_t *glob,
								fclaw2d_patch_t *patch,
								double alpha)
{
	/* We use the pack size here to make sure we are setting
	   everything correctly;  it isn't needed for memory
	   allocation */
	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int meqn = clawpatch_opt->meqn;
	int mbc = clawpatch_opt->mbc;
	int mint = clawpatch_opt->interp_stencil_width/2+1;  

	int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
	int wg = mx*my;  /* Whole grid but no ghost cells.  
						Ghost cells will be averaged from finer level. */
#if PATCH_DIM == 3
	int mz = clawpatch_opt->d3->mz;
	hole *= mz;
	wg *= mz;
#endif
	FCLAW_ASSERT(hole >= 0);

	int psize = (wg - hole)*meqn;
	FCLAW_ASSERT(psize > 0);

	/* Store time interpolated data that will be use in coarse grid
	   exchanges */
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	double *qlast = cp->griddata_last.dataPtr();
	double *qcurr = cp->griddata.dataPtr();
	double *qinterp = cp->griddata_time_interpolated.dataPtr();

	int ierror;

	/* Do interpolation only on interior, since ghost cells in qcurr
	   are invalid and will lead to floating point exceptions.
	   We do a ghost cell update at the intermediate time level.  The
	   neighboring fine grid will average to ghost cells of the interpolated
	   level, then the interpolated level is used to interpolate to fine grid
	   ghost cells. */
	fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
#if PATCH_DIM == 2
	clawpatch_vt->d2->fort_timeinterp(&mx,&my,&mbc,&meqn,&psize,
									  qcurr,qlast,qinterp,&alpha,&ierror);
#else
	clawpatch_vt->d3->fort_timeinterp(&mx,&my,&mz, &mbc,&meqn,&psize,
									  qcurr,qlast,qinterp,&alpha,&ierror);
#endif

}


/* ------------------------------------- Ghost filling  ------------------------------- */

static
void clawpatch_copy_face(fclaw2d_global_t *glob,
						 fclaw2d_patch_t *patch,
						 fclaw2d_patch_t *neighbor_patch,
						 int iface,
						 int time_interp,
						 fclaw2d_patch_transform_data_t *transform_data)

{

	const fclaw_clawpatch_options_t *clawpatch_opt = 
	            fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;

	/* This routine might be called between two time-sync patches */
	int meqn;
	double *qthis;
	fclaw2d_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);

	double *qneighbor;
	fclaw2d_clawpatch_timesync_data(glob,neighbor_patch,time_interp,&qneighbor,&meqn);

	if (fill_ghost(glob,time_interp))
	{
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
#if PATCH_DIM == 2		
		clawpatch_vt->d2->fort_copy_face(&mx,&my,&mbc,&meqn,qthis,
		                                 qneighbor,&iface,&transform_data);
#else
		int mz = clawpatch_opt->d3->mz;
		clawpatch_vt->d3->fort_copy_face(&mx,&my,&mz,&mbc,&meqn,qthis,
		                                 qneighbor,&iface,&transform_data);
#endif
	}
}

static
void clawpatch_average_face(fclaw2d_global_t *glob,
							fclaw2d_patch_t *coarse_patch,
							fclaw2d_patch_t *fine_patch,
							int idir,
							int iface_coarse,
							int p4est_refineFactor,
							int refratio,
							int time_interp,
							int igrid,
							fclaw2d_patch_transform_data_t* transform_data)
{
	int meqn;
	double *qcoarse;
	fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
	double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);


	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;

	const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
	int manifold = fclaw_opt->manifold;
	fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

	double *areacoarse = clawpatch_get_area(glob, coarse_patch);
	double *areafine = clawpatch_get_area(glob, fine_patch);

#if PATCH_DIM == 2
	/* These will be empty for non-manifolds cases */
	clawpatch_vt->d2->fort_average_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,
	                                	areacoarse,areafine, &idir,&iface_coarse, 
	                                	&p4est_refineFactor, &refratio,
										&igrid,&manifold,&transform_data);
#elif PATCH_DIM == 3
#if 0
	double *volcoarse = clawpatch_get_volume(glob, coarse_patch);
	double *volfine = clawpatch_get_volume(glob, fine_patch);
#endif
	/* These will be empty for non-manifolds cases */
	double *volcoarse = areacoarse;
	double *volfine = areafine;
	int mz = clawpatch_opt->d3->mz;
	clawpatch_vt->d3->fort_average_face(&mx,&my,&mz,&mbc,&meqn,qcoarse,qfine,
	                                    volcoarse,volfine, &idir,&iface_coarse, 
	                                    &p4est_refineFactor, &refratio,
	                                    &igrid,&manifold,&transform_data);
#endif
}

static
void clawpatch_interpolate_face(fclaw2d_global_t *glob,
								fclaw2d_patch_t *coarse_patch,
								fclaw2d_patch_t *fine_patch,
								int idir,
								int iface_coarse,
								int p4est_refineFactor,
								int refratio,
								int time_interp,
								int igrid,
								fclaw2d_patch_transform_data_t* transform_data)
{

	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

	int meqn;
	double *qcoarse;
	fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
	double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;

	if (fill_ghost(glob,time_interp))
	{
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
#if PATCH_DIM == 2
		clawpatch_vt->d2->fort_interpolate_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,&idir,
		                                        &iface_coarse, &p4est_refineFactor,
		                                        &refratio,&igrid,&transform_data);
#else
	int mz = clawpatch_opt->d3->mz;
	clawpatch_vt->d3->fort_interpolate_face(&mx,&my,&mz, &mbc,&meqn,qcoarse,qfine,
	                                        &idir,&iface_coarse, &p4est_refineFactor,
	                                        &refratio,&igrid,&transform_data);
#endif
	}
}

static
void clawpatch_copy_corner(fclaw2d_global_t *glob,
						   fclaw2d_patch_t *patch,
						   fclaw2d_patch_t *corner_patch,
						   int coarse_blockno,
						   int fine_blockno,
						   int icorner,
						   int time_interp,
						   fclaw2d_patch_transform_data_t *transform_data)
{

	const fclaw_clawpatch_options_t *clawpatch_opt = 
        	     fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;

	int meqn;
	double *qthis;
	fclaw2d_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);

	double *qcorner;
	fclaw2d_clawpatch_timesync_data(glob,corner_patch,time_interp,&qcorner,&meqn);

	if (fill_ghost(glob,time_interp))
	{
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
#if PATCH_DIM == 2
		clawpatch_vt->d2->fort_copy_corner(&mx,&my,&mbc,&meqn,qthis,qcorner,
		                                   &icorner, &transform_data);
#else
		int mz = clawpatch_opt->d3->mz;
		clawpatch_vt->d3->fort_copy_corner(&mx,&my,&mz, &mbc,&meqn,qthis,qcorner,
		                                   &icorner, &transform_data);
#endif
	}
}


static
void clawpatch_average_corner(fclaw2d_global_t *glob,
							  fclaw2d_patch_t *coarse_patch,
							  fclaw2d_patch_t *fine_patch,
							  int coarse_blockno,
							  int fine_blockno,
							  int coarse_corner,
							  int time_interp,
							  fclaw2d_patch_transform_data_t* transform_data)
{
	int meqn;
	double *qcoarse;
	fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

	double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	int manifold = fclaw_opt->manifold;
	if (fill_ghost(glob,time_interp))
	{
		int refratio = 2;
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

		double *areacoarse = clawpatch_get_area(glob, coarse_patch);
		double *areafine = clawpatch_get_area(glob, fine_patch);

#if PATCH_DIM == 2
   	    /* These will be empty for non-manifolds cases */
		clawpatch_vt->d2->fort_average_corner(&mx,&my,&mbc,&meqn,
		                                      &refratio,qcoarse,qfine,
		                                      areacoarse,areafine,
		                                      &manifold,&coarse_corner,&transform_data);
#elif PATCH_DIM == 3
     	/* These will be empty for non-manifolds cases */
#if 0     	
		double *volcoarse = clawpatch_get_volume(glob, coarse_patch);
		double *volfine = clawpatch_get_volume(glob, fine_patch);
#endif
		double *volcoarse = areacoarse;
		double *volfine = areafine;
		int mz = clawpatch_opt->d3->mz;		
		clawpatch_vt->d3->fort_average_corner(&mx,&my,&mz,&mbc,&meqn,
		                                      &refratio,qcoarse,qfine,
		                                      volcoarse,volfine,
		                                      &manifold,&coarse_corner,&transform_data);
#endif
	}
}

static
void clawpatch_interpolate_corner(fclaw2d_global_t* glob,
								  fclaw2d_patch_t* coarse_patch,
								  fclaw2d_patch_t* fine_patch,
								  int coarse_blockno,
								  int fine_blockno,
								  int coarse_corner,
								  int time_interp,
								  fclaw2d_patch_transform_data_t* transform_data)

{
	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;

	int meqn;
	double *qcoarse;
	fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

	double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

	if (fill_ghost(glob,time_interp))
	{
		int refratio = 2;
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
#if PATCH_DIM == 2
		clawpatch_vt->d2->fort_interpolate_corner(&mx,&my,&mbc,&meqn,
										          &refratio,qcoarse,qfine,
										          &coarse_corner,&transform_data);	
#else
		int mz = clawpatch_opt->d3->mz;
		clawpatch_vt->d3->fort_interpolate_corner(&mx,&my,&mz,&mbc,&meqn,
										          &refratio,qcoarse,qfine,
										          &coarse_corner,&transform_data);	
#endif

	}

}



/* -------------------------------- Regridding functions ------------------------------ */

static
int clawpatch_tag4refinement(fclaw2d_global_t *glob,
							 fclaw2d_patch_t *patch,
							 int blockno, int patchno,
							 int initflag)
{
	int meqn;
	double *q;
	fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	double refine_threshold = fclaw_opt->refine_threshold;

	int tag_patch;
	if (refine_threshold < 0) 
	{
		/* Always refine */
		tag_patch = 1;
	}
	else
	{
		tag_patch = 0;	

		/* This allows the user to specify a "exceeds_th" */
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

		int mx,my,mbc;
		double xlower,ylower,dx,dy;
		fclaw2d_global_set_global(glob);
#if PATCH_DIM == 2
		fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
		                            &xlower,&ylower,&dx,&dy);
		clawpatch_vt->d2->fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,
		                                      &dx,&dy, &blockno, q,
		                                      &refine_threshold,
		                                      &initflag,&tag_patch);
#elif PATCH_DIM == 3
		int mz;
		double zlower,dz;
		fclaw3dx_clawpatch_grid_data(glob,patch,&mx,&my,&mz, &mbc,
		                            &xlower,&ylower,&zlower, &dx,&dy,&dz);

		clawpatch_vt->d3->fort_tag4refinement(&mx,&my,&mz, &mbc,&meqn,
		                                      &xlower,&ylower,&zlower,
		                                      &dx,&dy, &dz, &blockno, q,
		                                      &refine_threshold,
		                                      &initflag,&tag_patch);
#endif		
		fclaw2d_global_unset_global();
	}
	return tag_patch;
}

static
int clawpatch_tag4coarsening(fclaw2d_global_t *glob,
							 fclaw2d_patch_t *fine_patches,
							 int blockno,
							 int patchno,
							 int initflag)
{
#if REFINE_DIM == 3
	fclaw_global_essentialf("tag4coarsening : Not implemented for full 3d " \
	                        "refinement.\n");
	exit(0);
#endif

	int mx,my,mbc,meqn;
	double xlower[4],ylower[4],dx,dy;
#if PATCH_DIM == 3
	int mz;
	double zlower, dz;
#endif

	double *q[4];  /* Only need four grids, even for extruded mesh case */
	for (int igrid = 0; igrid < 4; igrid++)
	{
		fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
#if PATCH_DIM == 2
		fclaw2d_clawpatch_grid_data(glob,&fine_patches[igrid],&mx,&my,&mbc,
		                            &xlower[igrid],&ylower[igrid],&dx,&dy);
#elif PATCH_DIM == 3
		/* For extruded meshes, zlower is the same for all patches. */
		fclaw2d_clawpatch_grid_data(glob,&fine_patches[igrid],&mx,&my,&mz,&mbc,
		                            &xlower[igrid],&ylower[igrid],&zlower,&dx,&dy,&dz);
#endif
	}

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	double coarsen_threshold = fclaw_opt->coarsen_threshold;

	int tag_patch = 0;
	if (coarsen_threshold > 0) 
	{		
		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
		fclaw2d_global_set_global(glob);
#if PATCH_DIM == 2
		/* Get xlower,ylower for each grid. */
		clawpatch_vt->d2->fort_tag4coarsening(&mx,&my,&mbc,&meqn,
		                                      xlower,ylower,&dx,&dy,
		                                      &blockno, q[0],q[1],q[2],q[3],
		                                      &coarsen_threshold,&initflag,&tag_patch);
#elif PATCH_DIM == 3

		clawpatch_vt->d3->fort_tag4coarsening(&mx,&my,&mz,&mbc,&meqn,
		                                      xlower,ylower,&zlower,&dx,&dy,&dz,
		                                      &blockno, q[0],q[1],q[2],q[3],
		                                      &coarsen_threshold,&initflag,&tag_patch);
#endif
		fclaw2d_global_unset_global();
	}
	else
	{
		/* Never coarsen */
	}
	return tag_patch == 1;
}

static
void clawpatch_interpolate2fine(fclaw2d_global_t* glob,
								fclaw2d_patch_t *coarse_patch,
								fclaw2d_patch_t* fine_patches,
								int this_blockno, int coarse_patchno,
								int fine0_patchno)
{

	const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;
	int meqn = clawpatch_opt->meqn;

	/* Loop over four siblings (z-ordering) */
	for (int igrid = 0; igrid < FCLAW2D_NUMSIBLINGS; igrid++)
	{
		fclaw2d_patch_t *fine_patch = &fine_patches[igrid];
		double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

		const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

		double *areacoarse = clawpatch_get_area(glob, coarse_patch);
		double *areafine = NULL;

		double *qcoarse = fclaw2d_clawpatch_get_q(glob,coarse_patch);

#if PATCH_DIM == 2
		if (fclaw_opt->manifold)
			areafine = clawpatch_get_area(glob, fine_patch);

		clawpatch_vt->d2->fort_interpolate2fine(&mx,&my,&mbc,&meqn,qcoarse,qfine,
											    areacoarse, areafine, &igrid,
											    &fclaw_opt->manifold);
#elif PATCH_DIM == 3
		double *volfine = areafine;
		if (fclaw_opt->manifold)
			 volfine = clawpatch_get_volume(glob, fine_patch);
		double *volcoarse = areacoarse;

		int mz = clawpatch_opt->d3->mz;
		clawpatch_vt->d3->fort_interpolate2fine(&mx,&my,&mz, &mbc,&meqn,qcoarse,qfine,
											    volcoarse, volfine, &igrid,
											    &fclaw_opt->manifold);
#endif
	}
}

static
void clawpatch_average2coarse(fclaw2d_global_t *glob,
							  fclaw2d_patch_t *fine_patches,
							  fclaw2d_patch_t *coarse_patch,
							  int blockno, int fine0_patchno,
							  int coarse_patchno)

{
	const fclaw_clawpatch_options_t *clawpatch_opt = 
	                fclaw_clawpatch_get_options(glob);	
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;
	int meqn = clawpatch_opt->meqn;

	double *areacoarse = clawpatch_get_area(glob, coarse_patch);

	double *qcoarse = fclaw2d_clawpatch_get_q(glob,coarse_patch);

	for(int igrid = 0; igrid < FCLAW2D_NUMSIBLINGS; igrid++)
	{
		fclaw2d_patch_t *fine_patch = &fine_patches[igrid];
		double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

		const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

		fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);


#if PATCH_DIM == 2	
		double *areafine =  NULL;
		if (fclaw_opt->manifold)
			areafine = clawpatch_get_area(glob, fine_patch);
		clawpatch_vt->d2->fort_average2coarse(&mx,&my,&mbc,&meqn,qcoarse,qfine,
										      areacoarse, areafine, &igrid,
										      &fclaw_opt->manifold);
#elif PATCH_DIM == 3
		double *volfine = NULL;
		double *volcoarse = areacoarse;
		if (fclaw_opt->manifold)
			volfine = clawpatch_get_volume(glob, fine_patch);
		int mz = clawpatch_opt->d3->mz;
		clawpatch_vt->d3->fort_average2coarse(&mx,&my,&mz, &mbc,&meqn,qcoarse,qfine,
										      volcoarse, volfine, &igrid,
										      &fclaw_opt->manifold);
#endif

	}
}

/* ------------------------------ Parallel ghost patches ------------------------------ */

/* This is called just to get a count of how much to pack */
static
size_t clawpatch_ghost_pack_elems(fclaw2d_global_t* glob)
{
	const fclaw_clawpatch_options_t *clawpatch_opt = 
					     	fclaw_clawpatch_get_options(glob);

	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;
	int meqn = clawpatch_opt->meqn;

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	int refratio = fclaw_opt->refratio;
	int packextra = fclaw_opt->ghost_patch_pack_numextrafields;
	int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
	int packregisters = fclaw_opt->time_sync;

	int mint = refratio*mbc;
	int nghost = mbc;

	/* Include size of conservation registers.  Save fluxes on each size, 
	   even though only one or two sides may be used. */
	// int frsize = 12*meqn*(mx + my); 

	int frsize = packregisters ? 2*(4*meqn+2)*(mx + my) : 0;
#if PATCH_DIM == 3
	int mz = clawpatch_opt->d3->mz;
	if (packregisters)
	{		
		fclaw_global_essentialf("clawpatch_ghost_comm: Conservation fix not yet " \
		                        "implemented in 3d\n");
		exit(0);
	}
#endif

	int wg = (2*nghost + mx)*(2*nghost + my);  /* Whole grid     */
	int hole = (mx - 2*mint)*(my - 2*mint);    /* Hole in center */

#if PATCH_DIM == 3	
	wg *= (mz + 2*nghost);
	hole *= (mz + 2*nghost); 
#endif
	FCLAW_ASSERT(hole >= 0);

	size_t psize = (wg - hole)*(meqn + packarea + packextra) + frsize;
	FCLAW_ASSERT(psize >= 0);

	return psize;
}    


static
void clawpatch_ghost_comm(fclaw2d_global_t* glob,
						  fclaw2d_patch_t* patch,
						  void *unpack_from_here, int time_interp,
						  int packmode)
{
	const fclaw_clawpatch_options_t *clawpatch_opt = 
	                        fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;
	int meqn = clawpatch_opt->meqn;

	int packarea = packmode/2;   // (0,1)/2 = 0;  (2,3)/2 = 1;

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	int packextra = fclaw_opt->ghost_patch_pack_numextrafields;
	int packregisters = fclaw_opt->time_sync;  // For conservation
	int refratio = fclaw_opt->refratio;

	int mint = refratio*mbc;   /* # interior cells needed for averaging */
	int nghost = mbc;          /* # ghost values needed for interpolation */

	/* Include size of conservation registers.  Save fluxes on each size, 
	   even though only one or two sides may be used.

	   y-face : gp,gm,f(qfront),f(qback) + edge-length + area : 4*meqn + 2
	   x-face : fp,fm,f(ql),f(qr) + edge-length + area        : 4*meqn + 2
	
	   Multiply everything by mz to get fluxes on each face of cube. 
	   (Time sync not yet implemented in 3d, though). 

	   */
	int frsize = packregisters ? 2*(4*meqn+2)*(mx + my) : 0;

#if PATCH_DIM == 3
	int mz = clawpatch_opt->d3->mz;
	if (packregisters)
	{		
		fclaw_global_essentialf("clawpatch_ghost_comm: Conservation fix not yet " \
		                        "implemented in 3d\n");
		exit(0);
	}
#endif

	/* wg   : whole grid
	   hole : part we don't store */ 
	/* This is computed twice - here, and in fclaw2d_clawpatch_ghost_packsize */
	int wg = (2*nghost + mx)*(2*nghost + my);
	int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */

#if PATCH_DIM == 3	
	wg *= (mz + 2*nghost);
	hole *= (mz + 2*nghost);  
#endif
	FCLAW_ASSERT(hole >= 0);

	size_t psize = (wg - hole)*(meqn + packarea + packextra) + frsize;
	FCLAW_ASSERT(psize > 0);

	size_t psize_check = clawpatch_ghost_pack_elems(glob);

	/* Check with routine that is used to allocate space */
	FCLAW_ASSERT(psize == psize_check);

	int ierror;
	fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

	double *qthis;
	fclaw2d_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);
	double *qpack = (double*) unpack_from_here;
#if PATCH_DIM == 2
	int qareasize = (wg - hole)*(meqn + packarea);
	double *area = clawpatch_get_area(glob, patch);	
	clawpatch_vt->d2->fort_local_ghost_pack(&mx,&my,&mbc,&meqn,&mint,qthis,area,
										    qpack,&qareasize,&packmode,&ierror);
#elif PATCH_DIM == 3
	int qvolsize = (wg - hole)*(meqn + packarea);
	double *vol = clawpatch_get_volume(glob, patch);	
	clawpatch_vt->d3->fort_local_ghost_pack(&mx,&my,&mz,&mbc,&meqn,&mint,qthis,vol,
										    qpack,&qvolsize,&packmode,&ierror);
#endif
	FCLAW_ASSERT(ierror == 0);

#if PATCH_DIM == 2
	qpack += qareasize;  /* Advance pointer */
#elif PATCH_DIM == 3
	qpack += qvolsize;  /* Advance pointer */
#endif

	int extrasize = (wg - hole)*(packextra);
	if (packextra)
	{
		FCLAW_ASSERT(extrasize > 0);
		FCLAW_ASSERT(clawpatch_vt->local_ghost_pack_aux != NULL);
		clawpatch_vt->local_ghost_pack_aux(glob,patch,mint,
		                                   qpack,extrasize,
		                                   packmode,&ierror);
		FCLAW_ASSERT(ierror == 0);
		qpack += extrasize; /* Advance pointer */
	}

	FCLAW_ASSERT(ierror == 0);
	if (packregisters)
	{
		FCLAW_ASSERT(frsize > 0);
		FCLAW_ASSERT(clawpatch_vt->time_sync_pack_registers != NULL);
		fclaw_clawpatch_packmode_t frpackmode = packmode % 2 == 0 ?  
		                                            CLAWPATCH_REGISTER_PACK : 
			                                        CLAWPATCH_REGISTER_UNPACK;
		clawpatch_vt->time_sync_pack_registers(glob, patch,
		                                       qpack,frsize,frpackmode,
		                                       &ierror);
		FCLAW_ASSERT(ierror == 0);
	}


	if (ierror > 0)
	{
		fclaw_global_essentialf("clawpatch_ghost_comm  : ierror = %d\n",ierror);
		exit(0);
	}
}


static size_t clawpatch_ghost_packsize(fclaw2d_global_t* glob)
{
	size_t esize = clawpatch_ghost_pack_elems(glob);
	return esize*sizeof(double);
}

static
void clawpatch_local_ghost_pack(fclaw2d_global_t *glob,
								fclaw2d_patch_t *patch,
								void *patch_data,
								int time_interp)
{
	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
	int packmode = 2*packarea;  // 0 or 2  (for pack)

	clawpatch_ghost_comm(glob,patch,patch_data, time_interp,packmode);
}

static
void clawpatch_remote_ghost_unpack(fclaw2d_global_t* glob,
								   fclaw2d_patch_t* patch,
								   int blockno,
								   int patchno,
								   void *qdata, int time_interp)
{
	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
	int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
	int packmode = 2*packarea + 1;  // 1 or 3  (for unpack)

	clawpatch_ghost_comm(glob,patch,qdata,time_interp,packmode);
}

static
void clawpatch_remote_ghost_build(fclaw2d_global_t *glob,
								  fclaw2d_patch_t *patch,
								  int blockno,
								  int patchno,
	                              fclaw2d_build_mode_t build_mode)
{
	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

	clawpatch_define(glob,patch,blockno,patchno,build_mode);

	if (fclaw_opt->manifold)
	{
		if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
		{
			/* Cell areas/volumes are not sent as MPI messages and so must
			   be recomputed
			*/
			//fclaw2d_metric_patch_build(glob,patch,blockno,patchno);
			fclaw2d_metric_patch_compute_area(glob,patch,blockno,patchno);
		}
	}
	/* Any metric terms we might need for the registers are packed */
#if 0	
    /* Setup for conservation correction */
	fclaw2d_clawpatch_time_sync_setup(glob,this_patch,blockno,patchno);
#endif	
}

static
void clawpatch_remote_ghost_delete(void *patchcp)
{
	clawpatch_delete(patchcp);
}

/* ---------------------------- Parallel partitioning --------------------------------- */

static
size_t clawpatch_partition_packsize(fclaw2d_global_t* glob)
{
	const fclaw_clawpatch_options_t *clawpatch_opt 
							  = fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int mbc = clawpatch_opt->mbc;
	int meqn = clawpatch_opt->meqn;
	size_t psize = meqn*(2*mbc + mx)*(2*mbc + my);  /* Store area */

#if PATCH_DIM == 3
	int mz = clawpatch_opt->d3->mz;
	psize *= (2*mbc + mz);
#endif

	return psize*sizeof(double);
}

static
void clawpatch_partition_pack(fclaw2d_global_t *glob,
							  fclaw2d_patch_t *patch,
							  int blockno,
							  int patchno,
							  void *pack_data_here)
	{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	FCLAW_ASSERT(cp != NULL);

	cp->griddata.copyToMemory((double*) pack_data_here);
}

static
void clawpatch_partition_unpack(fclaw2d_global_t *glob,  
								fclaw2d_domain_t *new_domain,
								fclaw2d_patch_t *patch,
								int blockno,
								int patchno,
								void *unpack_data_from_here)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);

	/* Time interp is false, since we only partition when all grids
	   are time synchronized and all flux registers are set to 
	   zero.  After copying data, we re-build patch with any 
	   data needed.  */
	cp->griddata.copyFromMemory((double*)unpack_data_from_here);
}

/* ------------------------------------ Virtual table  -------------------------------- */

static
fclaw_clawpatch_vtable_t* clawpatch_vt_new()
{
	fclaw_clawpatch_vtable_t * vt = FCLAW_ALLOC_ZERO (fclaw_clawpatch_vtable_t, 1);
#	if PATCH_DIM == 2
	vt->dim = 2;
	vt->d2 = FCLAW_ALLOC_ZERO (struct fclaw_clawpatch_vtable_2d, 1);
#	else
	vt->dim = 3;
	vt->d3 = FCLAW_ALLOC_ZERO (struct fclaw_clawpatch_vtable_3d, 1);
#	endif
    return vt;
}

static
void clawpatch_vt_destroy(void* vt)
{
#	if PATCH_DIM == 2
	FCLAW_FREE (((fclaw_clawpatch_vtable_t*) vt)->d2);
#	else
	FCLAW_FREE (((fclaw_clawpatch_vtable_t*) vt)->d3);
#	endif
    FCLAW_FREE (vt);
}

void fclaw2d_clawpatch_vtable_initialize(fclaw2d_global_t* glob, 
                                         int claw_version)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);

	fclaw2d_metric_vtable_initialize(glob);

	fclaw_clawpatch_vtable_t *clawpatch_vt = clawpatch_vt_new();

	/* Patch setup */
	patch_vt->patch_new             = clawpatch_new;
	patch_vt->patch_delete          = clawpatch_delete;
	patch_vt->build                 = clawpatch_build;
	patch_vt->build_from_fine       = clawpatch_build_from_fine;    

	/* Time stepping */
	patch_vt->restore_step          = clawpatch_restore_step;
	patch_vt->save_step             = clawpatch_save_step;
	patch_vt->setup_timeinterp      = clawpatch_setup_timeinterp;

	/* Ghost filling */
	patch_vt->copy_face            = clawpatch_copy_face;
	patch_vt->average_face         = clawpatch_average_face;
	patch_vt->interpolate_face     = clawpatch_interpolate_face;

	patch_vt->copy_corner          = clawpatch_copy_corner;
	patch_vt->average_corner       = clawpatch_average_corner;
	patch_vt->interpolate_corner   = clawpatch_interpolate_corner;

	/* Assume regular block corners;  Change by calling 'fclaw2d_clawpatch_use_pillowsphere' */    
	patch_vt->copy_block_corner          = clawpatch_copy_corner;
	patch_vt->average_block_corner       = clawpatch_average_corner;
	patch_vt->interpolate_block_corner   = clawpatch_interpolate_corner;

#if PATCH_DIM == 2
	/* Timing syncing module for conservation */
	patch_vt->time_sync_f2c       = fclaw2d_clawpatch_time_sync_f2c;
	patch_vt->time_sync_samesize  = fclaw2d_clawpatch_time_sync_samesize;
	patch_vt->time_sync_reset     = fclaw2d_clawpatch_time_sync_reset;
#endif

	/* Transform functions (defined in forestclaw2d */
	patch_vt->transform_init_data  = fclaw2d_clawpatch_transform_init_data;
	patch_vt->transform_face       = fclaw2d_clawpatch_face_transformation;       /* forestclaw2d.c */
	patch_vt->transform_face_intra = fclaw2d_clawpatch_face_transformation_intra; /* forestclaw2d.c */

	/* Regridding  functions */
	patch_vt->tag4refinement       = clawpatch_tag4refinement;
	patch_vt->tag4coarsening       = clawpatch_tag4coarsening;

	patch_vt->average2coarse       = clawpatch_average2coarse;
	patch_vt->interpolate2fine     = clawpatch_interpolate2fine;

	/* ghost patch */
	patch_vt->ghost_packsize       = clawpatch_ghost_packsize;
	patch_vt->local_ghost_pack     = clawpatch_local_ghost_pack;
	patch_vt->remote_ghost_build   = clawpatch_remote_ghost_build;
	patch_vt->remote_ghost_unpack  = clawpatch_remote_ghost_unpack;
	patch_vt->remote_ghost_delete  = clawpatch_remote_ghost_delete;


	/* partitioning */
	patch_vt->partition_packsize   = clawpatch_partition_packsize;
	patch_vt->partition_pack       = clawpatch_partition_pack;
	patch_vt->partition_unpack     = clawpatch_partition_unpack;

	/* output functions */
	clawpatch_vt->time_header_ascii  = fclaw_clawpatch_time_header_ascii;
	clawpatch_vt->cb_output_ascii    = cb_clawpatch_output_ascii; 

	/* Metric access */
	patch_vt->metric_patch         = clawpatch_get_metric_patch;

#if PATCH_DIM == 2
	/* Ghost pack for registers (doesn't depend on clawpack version) */
	clawpatch_vt->time_sync_pack_registers = fclaw2d_clawpatch_time_sync_pack_registers;
#endif

	/* Tagging functions.  The default uses option 'refinement_criteria'. */
#if PATCH_DIM == 2
	clawpatch_vt->d2->fort_user_exceeds_threshold = NULL;
#else
	clawpatch_vt->d3->fort_user_exceeds_threshold = NULL;
#endif

	/* Fortran functions that depend on data layout (version 4.6 or 5.0) */

#if PATCH_DIM == 2
	if (claw_version == 4)
	{
		/* Clawpatch settings functions */
		clawpatch_vt->d2->fort_average2coarse        = FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE;
		clawpatch_vt->d2->fort_interpolate2fine      = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE;

		clawpatch_vt->d2->fort_tag4refinement        = FCLAW2D_CLAWPATCH46_FORT_TAG4REFINEMENT;
		clawpatch_vt->d2->fort_tag4coarsening        = FCLAW2D_CLAWPATCH46_FORT_TAG4COARSENING;

		/* output functions */
		clawpatch_vt->fort_header_ascii              = FCLAW2D_CLAWPATCH46_FORT_HEADER_ASCII;
		clawpatch_vt->d2->fort_output_ascii          = FCLAW2D_CLAWPATCH46_FORT_OUTPUT_ASCII;

		/* Diagnostic functions */
		clawpatch_vt->conservation_check             = fclaw_clawpatch_diagnostics_cons_default;
		clawpatch_vt->compute_error                  = fclaw_clawpatch_diagnostics_error_default;
		clawpatch_vt->d2->fort_compute_patch_error   = NULL;   /* User defined */
		clawpatch_vt->d2->fort_compute_error_norm    = FCLAW2D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM;
		clawpatch_vt->d2->fort_compute_patch_area    = FCLAW2D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA;
		clawpatch_vt->d2->fort_conservation_check    = FCLAW2D_CLAWPATCH46_FORT_CONSERVATION_CHECK;

		/* Ghost cell exchange functions */
		clawpatch_vt->d2->fort_copy_face             = FCLAW2D_CLAWPATCH46_FORT_COPY_FACE;
		clawpatch_vt->d2->fort_average_face          = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE;
		clawpatch_vt->d2->fort_interpolate_face      = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE;

		clawpatch_vt->d2->fort_copy_corner           = FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER;
		clawpatch_vt->d2->fort_average_corner        = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER;
		clawpatch_vt->d2->fort_interpolate_corner    = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER;

		clawpatch_vt->local_ghost_pack_aux           = NULL;
		clawpatch_vt->d2->fort_local_ghost_pack      = FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;

		clawpatch_vt->d2->fort_timeinterp            = FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP;
	}
	else
	{
		/* Clawpatch settings functions */
		clawpatch_vt->d2->fort_average2coarse     = FCLAW2D_CLAWPATCH5_FORT_AVERAGE2COARSE;
		clawpatch_vt->d2->fort_interpolate2fine   = FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE2FINE;

		clawpatch_vt->d2->fort_tag4refinement     = FCLAW2D_CLAWPATCH5_FORT_TAG4REFINEMENT;
		clawpatch_vt->d2->fort_tag4coarsening     = FCLAW2D_CLAWPATCH5_FORT_TAG4COARSENING;

		/* output functions */
		clawpatch_vt->fort_header_ascii           = FCLAW2D_CLAWPATCH5_FORT_HEADER_ASCII;
		clawpatch_vt->d2->fort_output_ascii       = FCLAW2D_CLAWPATCH5_FORT_OUTPUT_ASCII;

		/* Diagnostic functions */
		clawpatch_vt->conservation_check          = fclaw_clawpatch_diagnostics_cons_default;
		clawpatch_vt->compute_error               = fclaw_clawpatch_diagnostics_error_default;
		clawpatch_vt->d2->fort_compute_patch_error  = NULL;   /* User defined */
		clawpatch_vt->d2->fort_compute_error_norm = FCLAW2D_CLAWPATCH5_FORT_COMPUTE_ERROR_NORM;
		clawpatch_vt->d2->fort_compute_patch_area = FCLAW2D_CLAWPATCH5_FORT_COMPUTE_PATCH_AREA;
		clawpatch_vt->d2->fort_conservation_check = FCLAW2D_CLAWPATCH5_FORT_CONSERVATION_CHECK;

		/* Ghost cell exchange functions */
		clawpatch_vt->d2->fort_copy_face          = FCLAW2D_CLAWPATCH5_FORT_COPY_FACE;
		clawpatch_vt->d2->fort_average_face       = FCLAW2D_CLAWPATCH5_FORT_AVERAGE_FACE;
		clawpatch_vt->d2->fort_interpolate_face   = FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_FACE;

		clawpatch_vt->d2->fort_copy_corner        = FCLAW2D_CLAWPATCH5_FORT_COPY_CORNER;
		clawpatch_vt->d2->fort_average_corner     = FCLAW2D_CLAWPATCH5_FORT_AVERAGE_CORNER;
		clawpatch_vt->d2->fort_interpolate_corner = FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_CORNER;

		clawpatch_vt->local_ghost_pack_aux        = NULL;
		clawpatch_vt->d2->fort_local_ghost_pack   = FCLAW2D_CLAWPATCH5_FORT_LOCAL_GHOST_PACK;

		clawpatch_vt->d2->fort_timeinterp         = FCLAW2D_CLAWPATCH5_FORT_TIMEINTERP;
	}
#elif PATCH_DIM == 3
	/* Signatures (defined in 'typedefs') for 3d Fortran routines are different 
	   those used in 2d routines above */
	if (claw_version == 4)
	{
		/* Clawpatch settings functions */
		clawpatch_vt->d3->fort_average2coarse        = FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE;
		clawpatch_vt->d3->fort_interpolate2fine      = FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE;

		clawpatch_vt->d3->fort_tag4refinement        = FCLAW3DX_CLAWPATCH46_FORT_TAG4REFINEMENT;
		clawpatch_vt->d3->fort_tag4coarsening        = FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING;

		/* output functions */
		clawpatch_vt->fort_header_ascii              = FCLAW3DX_CLAWPATCH46_FORT_HEADER_ASCII;
		clawpatch_vt->d3->fort_output_ascii          = FCLAW3DX_CLAWPATCH46_FORT_OUTPUT_ASCII;

		/* Ghost cell exchange functions */
		clawpatch_vt->d3->fort_copy_face             = FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE;
		clawpatch_vt->d3->fort_average_face          = FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE;
		clawpatch_vt->d3->fort_interpolate_face      = FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE;

		clawpatch_vt->d3->fort_copy_corner           = FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER;
		clawpatch_vt->d3->fort_average_corner        = FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER;
		clawpatch_vt->d3->fort_interpolate_corner    = FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER;

		clawpatch_vt->local_ghost_pack_aux           = NULL;
		clawpatch_vt->d3->fort_local_ghost_pack      = FCLAW3DX_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;

		clawpatch_vt->d3->fort_timeinterp            = FCLAW3DX_CLAWPATCH46_FORT_TIMEINTERP;

	}
	else if (claw_version == 5)
	{
		fclaw_global_essentialf("clawpatch_vtable_initialize : Version 5 not yet " \
		                        "implemented\n");
		exit(0);
	}	
#endif


	fclaw_clawpatch_diagnostics_vtable_initialize(glob);

	/* Set the virtual table, even if it isn't used */
#if PATCH_DIM == 2
	fclaw2d_clawpatch_pillow_vtable_initialize(glob, claw_version);
#else
	fclaw3d_clawpatch_pillow_vtable_initialize(glob, claw_version);
#endif

	clawpatch_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables, "fclaw_clawpatch") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fclaw_clawpatch", clawpatch_vt, clawpatch_vt_destroy);
}


/* ------------------------------- Public access functions ---------------------------- */

/* These functions are not virtualized and are not defined by the 
   patch interface */

#if PATCH_DIM == 2
fclaw_clawpatch_vtable_t* fclaw_clawpatch_vt(fclaw2d_global_t* glob)
{

	fclaw_clawpatch_vtable_t* clawpatch_vt = (fclaw_clawpatch_vtable_t*) 
	   											fclaw_pointer_map_get(glob->vtables, "fclaw_clawpatch");
	FCLAW_ASSERT(clawpatch_vt != nullptr);
	FCLAW_ASSERT(clawpatch_vt->is_set != 0);
	return clawpatch_vt;
}
#endif

/* Called from clawpack 4.6 and 5.0 */
void fclaw2d_clawpatch_save_current_step(fclaw2d_global_t* glob,
										 fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	cp->griddata_last = cp->griddata;
}


#if PATCH_DIM == 2
fclaw_clawpatch_t* 
fclaw2d_clawpatch_get_clawpatch(fclaw2d_patch_t* patch)
{
	return get_clawpatch(patch);
}
#endif


fclaw2d_metric_patch_t* 
fclaw2d_clawpatch_get_metric_patch(fclaw2d_patch_t* patch)
{
	return get_metric_patch(patch);
}

#if PATCH_DIM == 2
void fclaw2d_clawpatch_grid_data(fclaw2d_global_t* glob,
								 fclaw2d_patch_t* patch,
								 int* mx, int* my, int* mbc,
								 double* xlower, double* ylower,
								 double* dx, double* dy)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	*mx = cp->d2->mx;
	*my = cp->d2->my;
	*mbc = cp->mbc;
	*xlower = cp->d2->xlower;
	*ylower = cp->d2->ylower;
	*dx = cp->d2->dx;
	*dy = cp->d2->dy;
}
#elif PATCH_DIM == 3
void fclaw3dx_clawpatch_grid_data(fclaw2d_global_t* glob,
								 fclaw2d_patch_t* patch,
								 int* mx, int* my, int* mz, int* mbc,
								 double* xlower, double* ylower,
								 double* zlower, 
								 double* dx, double* dy, double* dz)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	*mx = cp->d3->mx;
	*my = cp->d3->my;
	*mz = cp->d3->mz;
	*mbc = cp->mbc;
	*xlower = cp->d3->xlower;
	*ylower = cp->d3->ylower;
	*zlower = cp->d3->zlower;
	*dx = cp->d3->dx;
	*dy = cp->d3->dy;
	*dz = cp->d3->dz;
}

/* The metric terms only know about fclaw3d routines;  not 3dx routines */
void fclaw3d_clawpatch_grid_data(fclaw2d_global_t* glob,
								 fclaw2d_patch_t* patch,
								 int* mx, int* my, int* mz, int* mbc,
								 double* xlower, double* ylower,
								 double* zlower, 
								 double* dx, double* dy, double* dz)

{
	fclaw3dx_clawpatch_grid_data(glob,patch,mx,my,mz,mbc,
	                             xlower,ylower,zlower,
	                             dx,dy,dz);
}

#endif


void fclaw2d_clawpatch_aux_data(fclaw2d_global_t *glob,
								fclaw2d_patch_t *patch,
								double **aux, int* maux)
{
	fclaw_clawpatch_t *cp = get_clawpatch (patch);

	*maux = cp->maux;
	*aux = cp->aux.dataPtr();
}

void fclaw2d_clawpatch_soln_data(fclaw2d_global_t* glob,
								 fclaw2d_patch_t* patch,
								 double **q, int* meqn)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	*q = cp->griddata.dataPtr();
	*meqn = cp->meqn;
}

void fclaw2d_clawpatch_rhs_data(fclaw2d_global_t* glob,
								 fclaw2d_patch_t* patch,
								 double **rhs, int *mfields)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	*rhs = cp->rhs.dataPtr();
	*mfields = cp->mfields;
}

void fclaw2d_clawpatch_elliptic_error_data(fclaw2d_global_t* glob,
                                           fclaw2d_patch_t* patch,
                                           double **err, int *mfields)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	*err = cp->elliptic_error.dataPtr();
	*mfields = cp->mfields;
}

void fclaw2d_clawpatch_elliptic_soln_data(fclaw2d_global_t* glob,
                                           fclaw2d_patch_t* patch,
                                           double **soln, int *mfields)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	*soln = cp->elliptic_soln.dataPtr();
	*mfields = cp->mfields;
}


double *fclaw2d_clawpatch_get_q(fclaw2d_global_t* glob,
								fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	return cp->griddata.dataPtr();
}

#if PATCH_DIM == 2
fclaw2d_clawpatch_registers_t* 
fclaw2d_clawpatch_get_registers(fclaw2d_global_t* glob,
                                  fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	return cp->d2->registers;
}
#endif


double* fclaw2d_clawpatch_get_error(fclaw2d_global_t* glob,
									fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	return cp->griderror.dataPtr();
}

double* fclaw2d_clawpatch_get_exactsoln(fclaw2d_global_t* glob,
									fclaw2d_patch_t* patch)
{
	fclaw_clawpatch_t *cp = get_clawpatch(patch);
	return cp->exactsolution.dataPtr();
}

void* fclaw2d_clawpatch_get_user_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->user_data;
}

void fclaw2d_clawpatch_set_user_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* patch,
                                    void *udata)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->user_data = udata;
}

void* fclaw2d_clawpatch_get_solver_data(fclaw2d_global_t* glob,
                                       fclaw2d_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->solver_data;
}

void fclaw2d_clawpatch_set_solver_data(fclaw2d_global_t* glob,
                                       fclaw2d_patch_t* patch,
                                       void *sdata)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->solver_data = sdata;
}

size_t fclaw2d_clawpatch_size(fclaw2d_global_t *glob)
{
	const fclaw_clawpatch_options_t *clawpatch_opt = 
					 fclaw_clawpatch_get_options(glob);
	int mx = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->mx : clawpatch_opt->d3->mx;
	int my = (clawpatch_opt->dim == 2) ? clawpatch_opt->d2->my : clawpatch_opt->d3->my;
	int meqn = clawpatch_opt->meqn;
	int mbc = clawpatch_opt->mbc;
	size_t size = (mx+2*mbc)*(my+2*mbc)*meqn;

#if PATCH_DIM == 3
	int mz = clawpatch_opt->d3->mz;
	size *= (mz + 2*mbc);
#endif

	return size;
}

void fclaw2d_clawpatch_timesync_data(fclaw2d_global_t* glob,
									 fclaw2d_patch_t* patch,
									 int time_interp,
									 double **q, int* meqn)
{
	fclaw_clawpatch_t* cp = get_clawpatch(patch);
	*q = q_time_sync(patch, time_interp);
	*meqn = cp->meqn;
}

double* fclaw2d_clawpatch_get_q_timesync(fclaw2d_global_t* glob,
										 fclaw2d_patch_t* patch,
										 int time_interp)
{
	return q_time_sync(patch, time_interp);
}

double* fclaw2d_clawpatch_get_area(fclaw2d_global_t* glob,
								   fclaw2d_patch_t* patch)
{
	return clawpatch_get_area(glob, patch);
}

#if PATCH_DIM == 2
void fclaw2d_clawpatch_metric_scalar(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* patch,
                                     double **area, double** edgelengths,
                                     double **curvature)
{
	fclaw2d_metric_patch_scalar(glob,patch,area,edgelengths,
	                            curvature);
}

void fclaw2d_clawpatch_metric_vector(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch,
                                     double **xnormals, double **ynormals,
                                     double **xtangents, double **ytangents,
                                     double **surfnormals)
{
	fclaw2d_metric_patch_vector(glob,patch,xnormals,ynormals,
	                            xtangents,ytangents,surfnormals);
}




void fclaw2d_clawpatch_metric_data(fclaw2d_global_t* glob,
								   fclaw2d_patch_t* patch,
								   double **xp, double **yp, double **zp,
								   double **xd, double **yd, double **zd,
								   double **area)
{
	fclaw2d_metric_patch_mesh_data(glob,patch,xp,yp,zp,xd,yd,zd,area);
}

void fclaw2d_clawpatch_metric_data2(fclaw2d_global_t* glob,
									fclaw2d_patch_t* patch,
									double **xnormals, double **ynormals,
									double **xtangents, double **ytangents,
									double **surfnormals,
									double **edgelengths, double **curvature)
{
	fclaw2d_metric_patch_mesh_data2(glob,patch,xnormals,ynormals,
									xtangents,ytangents,surfnormals,
									edgelengths,curvature);
}
#elif PATCH_DIM == 3

void fclaw3d_clawpatch_metric_scalar(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* patch,
                                     double **volume, double** faceareas)
{
	fclaw3d_metric_patch_scalar(glob,patch,volume,faceareas);
}

void fclaw2d_clawpatch_metric_basis(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch,
                                     double **xrot, double** yrot, double** zrot)
{
	fclaw3d_metric_patch_basis(glob,patch,xrot, yrot, zrot);
}




void fclaw3d_clawpatch_mesh_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* patch,
                                 double **xp, double **yp, double **zp,
                                 double **xd, double **yd, double **zd,
                                 double **volume, double **faceareas)
{
	fclaw3d_metric_patch_mesh_data(glob,patch,xp,yp,zp,xd,yd,zd,
	                               volume, faceareas);
}


#endif

