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

#include "fc3d_clawpack46.h"
#include "fc3d_clawpack46_options.h"
#include "fc3d_clawpack46_fort.h"

#include <fclaw_pointer_map.h>

#include <fclaw_clawpatch.hpp>
#include <fclaw_clawpatch.h>

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw3d_clawpatch_fort.h>

#ifdef FCLAW_ENABLE_HDF5
#include <fclaw_clawpatch_output_hdf5.h>
#endif

#include <fclaw3d_metric.h>
#include <fclaw3d_metric.hpp>


#include <fclaw_patch.h>
#include <fclaw_global.h>
#include <fclaw_vtable.h>
#include <fclaw_options.h>
#include <fclaw2d_defs.h>


/* --------------------- Clawpack solver functions (required) ------------------------- */

static
void clawpack46_setprob(fclaw_global_t *glob)
{
	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);
	if (claw46_vt->fort_setprob != NULL)
	{
		claw46_vt->fort_setprob();
	}
}


static
void clawpack46_qinit(fclaw_global_t *glob,
					  fclaw_patch_t *patch,
					  int blockno,
					  int patchno)
{
	int mx,my,mz,mbc;
	double dx,dy,dz, xlower,ylower, zlower;
	fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz,&mbc,
								&xlower,&ylower,&zlower, 
								&dx,&dy,&dz);

	int meqn;
	double *q;
	fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);

	FCLAW_ASSERT(claw46_vt->fort_qinit != NULL); /* Must be initialized */

	/* Call to classic Clawpack 'qinit' routine.  This must be user defined */
	FC3D_CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_qinit(&meqn,&mbc,&mx,&my,&mz, &xlower,&ylower,&zlower, 
	                      &dx,&dy,&dz, q, &maux,aux);
	FC3D_CLAWPACK46_UNSET_BLOCK();
}


static
void clawpack46_bc3(fclaw_global_t *glob,
					fclaw_patch_t *patch,
					int blockno,
					int patchno,
					double t,
					double dt,
					int intersects_phys_bdry[],
					int time_interp)
{
	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);

	fc3d_clawpack46_options_t *clawpack_options = fc3d_clawpack46_get_options(glob);

	FCLAW_ASSERT(claw46_vt->fort_bc3 != NULL);

	int mx,my,mz, mbc;
	double dx,dy,dz, xlower,ylower, zlower;
	fclaw_clawpatch_3d_grid_data(glob,patch, &mx,&my,&mz,&mbc,
								&xlower,&ylower,&zlower, &dx,&dy,&dz);

	double *aux;
	int maux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	int *block_mthbc = clawpack_options->mthbc;

	/* Set a local copy of mthbc that can be used for a patch. */
	int mthbc[6];
	if(glob->domain->refine_dim == 2)
	{
		for(int i = 0; i < 6; i++)
		{
			if (i < 4)			
				if (intersects_phys_bdry[i])
					mthbc[i] = block_mthbc[i];
				else
					mthbc[i] = -1;
			else
				mthbc[i] = block_mthbc[i];
		}
	}
	else 
	{
		for(int i = 0; i < 6; i++)
		{
			if (intersects_phys_bdry[i])
				mthbc[i] = block_mthbc[i];
			else
				mthbc[i] = -1;
		}
	}

	/*
	  We may be imposing boundary conditions on time-interpolated data;
	  and is being done just to help with fine grid interpolation.
	  In this case, this boundary condition won't be used to update
	  anything
	*/
	double *q;
	int meqn;
	fclaw_clawpatch_timesync_data(glob,patch,time_interp,&q,&meqn);

	FC3D_CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_bc3(&meqn,&mbc,&mx,&my,&mz,&xlower,&ylower,&zlower,
						&dx,&dy,&dz, q,&maux,aux,&t,&dt,mthbc);
	FC3D_CLAWPACK46_UNSET_BLOCK();
}


static
void clawpack46_b4step3(fclaw_global_t *glob,
						fclaw_patch_t *patch,
						int blockno,
						int patchno,
						double t, double dt)

{
	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);
	if (claw46_vt->fort_b4step3 == NULL)
		return;

	int mx,my,mz,mbc;
	double xlower,ylower,zlower,dx,dy,dz;
	fclaw_clawpatch_3d_grid_data(glob,patch, &mx,&my,&mz, &mbc,
								&xlower,&ylower,&zlower,&dx,&dy,&dz);

	int meqn;
	double *q;
	fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	FC3D_CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_b4step3(&mbc,&mx,&my,&mz, &meqn,q, 
	                        &xlower,&ylower,&zlower,
							&dx,&dy,&dz, &t,&dt,&maux,aux);
	FC3D_CLAWPACK46_UNSET_BLOCK();
}

static
void clawpack46_src3(fclaw_global_t *glob,
					 fclaw_patch_t *patch,
					 int blockno,
					 int patchno,
					 double t,
					 double dt)
{
	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);

	if (claw46_vt->fort_src3 == NULL)
		return;


	int mx,my,mz, mbc;
	double xlower,ylower,zlower,dx,dy,dz;
	fclaw_clawpatch_3d_grid_data(glob,patch, &mx,&my,&mz, &mbc,
								&xlower,&ylower,&zlower, &dx,&dy, &dz);

	double *q;
	int meqn;
	fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

	double *aux;
	int maux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	FC3D_CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_src3(&meqn,&mbc,&mx,&my,&mz, &xlower,&ylower,&zlower,
						 &dx,&dy,&dz, q,&maux,aux,&t,&dt);
	FC3D_CLAWPACK46_UNSET_BLOCK();
}


/* This can be used as a value for patch_vt->patch_setup */
static
void clawpack46_setaux(fclaw_global_t *glob,
					   fclaw_patch_t *patch,
					   int blockno,
					   int patchno)
	{
	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);

	if (claw46_vt->fort_setaux == NULL)
	{
		/* The user has not specified any aux routine */
		return;
	}

	if (fclaw_patch_is_ghost(patch))
	{
		/* This is going to be removed at some point */
		return;
	}

	int mx,my,mz,mbc;
	double xlower,ylower,zlower, dx,dy, dz;
	fclaw_clawpatch_3d_grid_data(glob,patch, &mx,&my,&mz,&mbc,
								&xlower,&ylower,&zlower, &dx,&dy, &dz);
	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	FC3D_CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_setaux(&mbc,&mx,&my,&mz, &xlower,&ylower,&zlower, 
	                       &dx,&dy,&dz, &maux,aux);
	FC3D_CLAWPACK46_UNSET_BLOCK();
}

/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
static
double clawpack46_step3(fclaw_global_t *glob,
						fclaw_patch_t *patch,
						int blockno,
						int patchno,
						double t,
						double dt)
{
	// const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

	fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);
	FCLAW_ASSERT(claw46_vt->fort_rpn3 != NULL);

	const fc3d_clawpack46_options_t* clawpack_options = 
	                fc3d_clawpack46_get_options(glob);

	if (clawpack_options->order[1] > 0)
		FCLAW_ASSERT(claw46_vt->fort_rpt3 != NULL);

	if (clawpack_options->order[2] > 0)
		FCLAW_ASSERT(claw46_vt->fort_rptt3 != NULL);

	int level = patch->level;

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	fclaw_clawpatch_save_current_step(glob, patch);

	int mx, my, mz, mbc;
	double xlower, ylower, zlower, dx,dy, dz;
	fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz,&mbc,
								&xlower,&ylower,&zlower, &dx,&dy,&dz);

	int meqn;
	double *qold;
	fclaw_clawpatch_soln_data(glob,patch,&qold,&meqn);


	int mwaves = clawpack_options->mwaves;

	int maxm = fmax(mx,fmax(my,mz));

	double cflgrid = 0.0;

#if 0		  
	fclaw3dx_clawpatch_registers_t* cr = 
		  fclaw3dx_clawpatch_get_registers(glob,patch);

	/* Evaluate fluxes needed in correction terms */
	if (fclaw_opt->time_sync && fclaw_opt->flux_correction)
	{
		FCLAW_ASSERT(claw46_vt->fort_rpn2_cons != NULL);
		double *qvec   = FCLAW_ALLOC(double, meqn);
		double *auxvec_center = FCLAW_ALLOC(double, maux);
		double *auxvec_edge = FCLAW_ALLOC(double, maux);
		double *flux   = FCLAW_ALLOC(double, meqn);     /* f(qr) - f(ql) = amdq+apdq */

		CLAWPACK46_TIME_SYNC_STORE_FLUX(&mx,&my,&mbc,&meqn,&maux,
		                                &blockno,&patchno, &dt,
										  cr->edgelengths[0], 
										  cr->edgelengths[1], 
										  cr->edgelengths[2], 
										  cr->edgelengths[3],
										  qold,aux,
										  cr->edge_fluxes[0],cr->edge_fluxes[1],
										  cr->edge_fluxes[2],cr->edge_fluxes[3],
										  claw46_vt->fort_rpn2_cons,
										  qvec,auxvec_center,auxvec_edge,flux);

		FCLAW_FREE(qvec);
		FCLAW_FREE(auxvec_center);
		FCLAW_FREE(auxvec_edge);
		FCLAW_FREE(flux);
	}
#endif	


#if 0
  	 /* from stepgrid */
    msize = (maxm + 2*mbc)
    mwork = msize*(46*meqn + (meqn+1)*mwaves + 9*maux + 3)
#endif	


	int msize = maxm + 2*mbc;
	int mwork = msize*(46*meqn + (meqn+1)*mwaves + 9*maux + 3);
	double* work = FCLAW_ALLOC(double,mwork);

	int size = meqn*(mx+2*mbc)*(my+2*mbc)*(mz + 2*mbc);
	double* fp = FCLAW_ALLOC(double,size);
	double* fm = FCLAW_ALLOC(double,size);
	double* gp = FCLAW_ALLOC(double,size);
	double* gm = FCLAW_ALLOC(double,size);
	double* hp = FCLAW_ALLOC(double,size);
	double* hm = FCLAW_ALLOC(double,size);

	int ierror = 0;
	int* block_corner_count = fclaw_patch_block_corner_count(glob,patch);

#if 0
	if (claw46_vt->flux2 == NULL)
	{
		claw46_vt->flux2 = (clawpack_options->use_fwaves != 0) ? &CLAWPACK46_FLUX2FW : 
		                       &CLAWPACK46_FLUX2;	
	}
#endif

	/* NOTE: qold will be overwritten in this step */
	FC3D_CLAWPACK46_SET_BLOCK(&blockno);
	CLAWPACK46_STEP3_WRAP(&maxm, &meqn, &maux, &mbc, clawpack_options->method,
						  clawpack_options->mthlim, &clawpack_options->mcapa,
						  &mwaves,&mx, &my, &mz, qold, aux, &dx, &dy, &dz, 
						  &dt, &cflgrid, work, &mwork, &xlower, &ylower, &zlower,
						  &level,&t, fp, fm, gp, gm, hp, hm, 
						  claw46_vt->fort_rpn3, claw46_vt->fort_rpt3,
						  claw46_vt->fort_rptt3,
						  &clawpack_options->use_fwaves, block_corner_count, 
						  &ierror);
	FC3D_CLAWPACK46_UNSET_BLOCK();

	FCLAW_ASSERT(ierror == 0);

#if 0
	if (fclaw_opt->time_sync && fclaw_opt->fluctuation_correction)
	{

		CLAWPACK46_TIME_SYNC_ACCUMULATE_WAVES(&mx,&my,&mbc,&meqn, &dt, &dx, 
		                                      &dy, &patchno,
		                                      cr->edgelengths[0],
		                                      cr->edgelengths[1],
		                                      cr->edgelengths[2],
		                                      cr->edgelengths[3],
		                                      fp,fm,gp,gm,
		                                      cr->fp[0],cr->fp[1],
		                                      cr->fm[0],cr->fm[1],
		                                      cr->gp[0],cr->gp[1],
		                                      cr->gm[0],cr->gm[1]);
	}
#endif			


	FCLAW_FREE(fp);
	FCLAW_FREE(fm);
	FCLAW_FREE(gp);
	FCLAW_FREE(gm);
	FCLAW_FREE(hp);
	FCLAW_FREE(hm);
	FCLAW_FREE(work);

	return cflgrid;
}

static
double clawpack46_update(fclaw_global_t *glob,
                         fclaw_patch_t *patch,
                         int blockno,
                         int patchno,
                         double t,
                         double dt, 
                         void* user)
{

    fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);
    if (claw46_vt->b4step3 != NULL)
    {
        fclaw_timer_start_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_B4STEP2]);               
        claw46_vt->b4step3(glob,
                           patch,
                           blockno,
                           patchno,t,dt);

        fclaw_timer_stop_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_B4STEP2]);               
    }

    fclaw_timer_start_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_STEP2]);       

    double maxcfl = clawpack46_step3(glob,
                                     patch,
                                     blockno,
                                     patchno,t,dt);

    fclaw_timer_stop_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_STEP2]);       

    const fc3d_clawpack46_options_t* clawpack_options = fc3d_clawpack46_get_options(glob);
    if (clawpack_options->src_term > 0 && claw46_vt->src3 != NULL)
    {
        claw46_vt->src3(glob,
                        patch,
                        blockno,
                        patchno,t,dt);
    }
    return maxcfl;
}


/* ---------------------------------- Output functions -------------------------------- */

static
void clawpack46_output(fclaw_global_t *glob, int iframe)
{
	const fc3d_clawpack46_options_t* clawpack_options 
	                  = fc3d_clawpack46_get_options(glob);
	if (clawpack_options->ascii_out != 0)
		fclaw_clawpatch_output_ascii(glob,iframe);

	if (clawpack_options->vtk_out != 0)
		fclaw_clawpatch_output_vtk(glob,iframe);

#ifdef FCLAW_ENABLE_HDF5
	if (clawpack_options->hdf_out != 0)
		fclaw_clawpatch_output_hdf5(glob,iframe);
#endif
}



/* ---------------------------------- Virtual table  ---------------------------------- */

static
fc3d_clawpack46_vtable_t* clawpack46_vt_new()
{
    return (fc3d_clawpack46_vtable_t*) FCLAW_ALLOC_ZERO (fc3d_clawpack46_vtable_t, 1);
}

static
void clawpack46_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

void fc3d_clawpack46_solver_initialize(fclaw_global_t* glob)
{
	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
	fc3d_clawpack46_options_t* clawopt = fc3d_clawpack46_get_options(glob);

	if(clawpatch_opt->patch_dim != 3)
	{
		fclaw_abortf("Clawpatch dimension set to 2d. fc3d_clawpack46 is only for 3d.");
	}

    clawopt->method[6] = clawpatch_opt->maux;

    if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack : bad maux/mcapa combination\n");
        exit(FCLAW_EXIT_ERROR);
    }

	int claw_version = 4;
	fclaw_clawpatch_vtable_initialize(glob, claw_version);
    //fclaw_clawpatch_vtable_t*      clawpatch_vt = fclaw_clawpatch_vt();

	fclaw_vtable_t*                fc_vt = fclaw_vt(glob);
	fclaw_patch_vtable_t*          patch_vt = fclaw_patch_vt(glob);  

	fc3d_clawpack46_vtable_t*  claw46_vt = clawpack46_vt_new();

	/* ForestClaw vtable items */
	fc_vt->output_frame      = clawpack46_output;
	fc_vt->problem_setup     = clawpack46_setprob;    

	/* These could be over-written by user specific settings */
	patch_vt->initialize                     = clawpack46_qinit;
	patch_vt->setup                          = clawpack46_setaux;  
	patch_vt->physical_bc                    = clawpack46_bc3;
	patch_vt->single_step_update             = clawpack46_update;

#if 0
	/* Conservation updates (based on Clawpack updates) */
	clawpatch_vt->fort_time_sync_f2c         = CLAWPACK46_FORT_TIME_SYNC_F2C;
	clawpatch_vt->fort_time_sync_samesize    = CLAWPACK46_FORT_TIME_SYNC_SAMESIZE;
#endif	

	/* Wrappers so that user can change argument list */
	claw46_vt->b4step3                       = clawpack46_b4step3;
	claw46_vt->src3                          = clawpack46_src3;

	/* Required functions  - error if NULL */
	claw46_vt->fort_bc3       = CLAWPACK46_BC3_DEFAULT;
	claw46_vt->fort_qinit     = NULL;
	claw46_vt->fort_rpn3      = NULL;
	claw46_vt->fort_rpt3      = NULL;
	claw46_vt->fort_rptt3     = NULL;

	/* Optional functions - call only if non-NULL */
	claw46_vt->fort_setprob   = NULL;
	claw46_vt->fort_setaux    = NULL;
	claw46_vt->fort_b4step3   = NULL;
	claw46_vt->fort_src3      = NULL;

	claw46_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fc3d_clawpack46") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fc3d_clawpack46", claw46_vt, clawpack46_vt_destroy);
}


/* ----------------------------- User access to solver functions --------------------------- */


/* These are here in case the user wants to call Clawpack routines directly */

fc3d_clawpack46_vtable_t* fc3d_clawpack46_vt(fclaw_global_t* glob)
{
	fc3d_clawpack46_vtable_t* claw46_vt = (fc3d_clawpack46_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fc3d_clawpack46");
	FCLAW_ASSERT(claw46_vt != NULL);
	FCLAW_ASSERT(claw46_vt->is_set != 0);
	return claw46_vt;
}

/* This should only be called when a new fclaw3dx_clawpatch_t is created. */
void fc3d_clawpack46_set_capacity(fclaw_global_t *glob,
								  fclaw_patch_t *patch,
								  int blockno,
								  int patchno)
{
	fc3d_clawpack46_options_t *clawopt = fc3d_clawpack46_get_options(glob);
	int mcapa = clawopt->mcapa;

	int mx,my,mz, mbc;
	double dx,dy,dz, xlower,ylower,zlower;
	fclaw_clawpatch_3d_grid_data(glob,patch, &mx,&my,&mz,&mbc,
								&xlower,&ylower,&zlower,&dx,&dy,&dz);

	double *volume = fclaw_clawpatch_get_3d_volume(glob,patch);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);
	FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

	FC3D_CLAWPACK46_SET_CAPACITY(&mx,&my,&mz,&mbc,&dx,&dy,&dz,volume,&mcapa,
	                             &maux,aux);
}



/* -------------------------- Public interface to Clawpack wrappers --------------------*/

/* These are overkill;  it isn't obvious why the user would want these */

void fc3d_clawpack46_setprob(fclaw_global_t *glob)
{
	clawpack46_setprob(glob);
}

/* This can be set to claw46_vt->src3 */
void fc3d_clawpack46_src3(fclaw_global_t* glob,
						  fclaw_patch_t *patch,
						  int blockno,
						  int patchno,
						  double t,
						  double dt)
{
	clawpack46_src3(glob,patch,blockno,blockno,t,dt);
}


void fc3d_clawpack46_setaux(fclaw_global_t *glob,
							fclaw_patch_t *patch,
							int blockno,
							int patchno)
{
	clawpack46_setaux(glob,patch,blockno,patchno);
}


void fc3d_clawpack46_qinit(fclaw_global_t *glob,
						   fclaw_patch_t *patch,
						   int blockno,
						   int patchno)
{
	clawpack46_qinit(glob,patch,blockno,patchno);
}

void fc3d_clawpack46_b4step3(fclaw_global_t* glob,
							 fclaw_patch_t *patch,
							 int blockno,
							 int patchno,
							 double t,
							 double dt)
{
	clawpack46_b4step3(glob,patch,blockno,patchno,t,dt);
}

void fc3d_clawpack46_bc3(fclaw_global_t *glob,
						 fclaw_patch_t *patch,
						 int blockno,
						 int patchno,
						 double t,
						 double dt,
						 int intersects_bc[],
						 int time_interp)
{
	clawpack46_bc3(glob,patch,blockno,blockno,t,dt,
				   intersects_bc,time_interp);
}














