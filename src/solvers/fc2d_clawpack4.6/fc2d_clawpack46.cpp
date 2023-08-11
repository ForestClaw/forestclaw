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

#include "fc2d_clawpack46.h"
#include "fc2d_clawpack46_options.h"
#include "fc2d_clawpack46_fort.h"

#include <fclaw_clawpatch.hpp>
#include <fclaw_clawpatch.h>

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fclaw2d_clawpatch_conservation.h>

#include <fclaw2d_patch.h>
#include <fclaw_global.h>
#include <fclaw_vtable.h>
#include <fclaw_options.h>
#include <fclaw2d_defs.h>

#include <fclaw_pointer_map.h>


/* --------------------- Clawpack solver functions (required) ------------------------- */

static
void clawpack46_setprob(fclaw_global_t *glob)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);
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
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);
	FCLAW_ASSERT(claw46_vt->fort_qinit != NULL); /* Must be initialized */

	int mx,my,mbc;
	double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);
	int meqn;
	double *q;
	fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	int maxmx = mx;
	int maxmy = my;

	/* Call to classic Clawpack 'qinit' routine.  This must be user defined */
	CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_qinit(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
						  &maux,aux);
	CLAWPACK46_UNSET_BLOCK();
}


static
void clawpack46_bc2(fclaw_global_t *glob,
					fclaw_patch_t *patch,
					int blockno,
					int patchno,
					double t,
					double dt,
					int intersects_phys_bdry[],
					int time_interp)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);
	FCLAW_ASSERT(claw46_vt->fort_bc2 != NULL);

	int mx,my,mbc;
	double xlower,ylower,dx,dy;
	fclaw2d_clawpatch_grid_data(glob,patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	int maux;
	double* aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	int maxmx = mx;
	int maxmy = my;

	fc2d_clawpack46_options_t *clawpack_options = fc2d_clawpack46_get_options(glob);
	int *block_mthbc = clawpack_options->mthbc;

	/* Set a local copy of mthbc that can be used for a patch. */
	int mthbc[4];
	for(int i = 0; i < 4; i++)
	{
		if (intersects_phys_bdry[i])
			mthbc[i] = block_mthbc[i];
		else
			mthbc[i] = -1;
	}

	/*
	  We may be imposing boundary conditions on time-interpolated data;
	  and is being done just to help with fine grid interpolation.
	  In this case, this boundary condition won't be used to update
	  anything
	*/
	int meqn;
	double *q;
	fclaw_clawpatch_timesync_data(glob,patch,time_interp,&q,&meqn);

	CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_bc2(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,
						&dx,&dy,q,&maux,aux,&t,&dt,mthbc);
	CLAWPACK46_UNSET_BLOCK();
}


static
void clawpack46_b4step2(fclaw_global_t *glob,
						fclaw_patch_t *patch,
						int blockno,
						int patchno,
						double t, double dt)

{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);
	if (claw46_vt->fort_b4step2 == NULL)
		return;

	int mx,my,mbc;
	double xlower,ylower,dx,dy;
	fclaw2d_clawpatch_grid_data(glob,patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	int meqn;
	double *q;
	fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	int maxmx = mx;
	int maxmy = my;

	CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_b4step2(&maxmx,&maxmy,&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
							&dx,&dy,&t,&dt,&maux,aux);
	CLAWPACK46_UNSET_BLOCK();
}

static
void clawpack46_src2(fclaw_global_t *glob,
					 fclaw_patch_t *patch,
					 int blockno,
					 int patchno,
					 double t,
					 double dt)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);

	if (claw46_vt->fort_src2 == NULL)
	{
		/* User has not set a fortran routine */
		return;
	}

	int mx,my,mbc;
	double xlower,ylower,dx,dy;
	fclaw2d_clawpatch_grid_data(glob,patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	int meqn;
	double *q;
	fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	int maxmx = mx;
	int maxmy = my;

	CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_src2(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,
						 &dx,&dy,q,&maux,aux,&t,&dt);
	CLAWPACK46_UNSET_BLOCK();
}


/* This can be used as a value for patch_vt->patch_setup */
static
void clawpack46_setaux(fclaw_global_t *glob,
					   fclaw_patch_t *patch,
					   int blockno,
					   int patchno)
	{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);

	if (claw46_vt->fort_setaux == NULL)
		return;

	if (fclaw2d_patch_is_ghost(patch))
	{
		/* This is going to be removed at some point */
		return;
	}

	int mx,my,mbc;
	double xlower,ylower,dx,dy;
	fclaw2d_clawpatch_grid_data(glob,patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	int maxmx = mx;
	int maxmy = my;

	CLAWPACK46_SET_BLOCK(&blockno);
	claw46_vt->fort_setaux(&maxmx,&maxmy,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
						   &maux,aux);
	CLAWPACK46_UNSET_BLOCK();
}

/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
static
double clawpack46_step2(fclaw_global_t *glob,
						fclaw_patch_t *patch,
						int blockno,
						int patchno,
						double t,
						double dt)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);
	const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
	const fc2d_clawpack46_options_t* clawpack_options;

	clawpack_options = fc2d_clawpack46_get_options(glob);


	if (clawpack_options->use_fwaves)
	{
		FCLAW_ASSERT(claw46_vt->fort_rpn2fw != NULL);
		if (clawpack_options->order[1] > 0)
			FCLAW_ASSERT(claw46_vt->fort_rpt2fw != NULL);
	}
	else
	{
		FCLAW_ASSERT(claw46_vt->fort_rpn2 != NULL);
		if (clawpack_options->order[1] > 0)
			FCLAW_ASSERT(claw46_vt->fort_rpt2 != NULL);
	}



	int level = patch->level;

	int maux;
	double *aux; 
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	fclaw_clawpatch_save_current_step(glob, patch);

	int mx, my, mbc;
	double xlower, ylower, dx,dy;
	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	int meqn;
	double *qold;
	fclaw_clawpatch_soln_data(glob,patch,&qold,&meqn);

	int mwaves = clawpack_options->mwaves;
	int maxm = fmax(mx,my);

	double cflgrid = 0.0;

	fclaw2d_clawpatch_registers_t* cr = 
		  fclaw2d_clawpatch_get_registers(glob,patch);

	int* block_corner_count = fclaw2d_patch_block_corner_count(glob,patch);

	/* Evaluate fluxes needed in correction terms */
	if (fclaw_opt->time_sync && fclaw_opt->flux_correction)
	{
		FCLAW_ASSERT(claw46_vt->fort_rpn2_cons != NULL);
		double *qvec          = FCLAW_ALLOC(double, meqn);
		double *auxvec_center = FCLAW_ALLOC(double, maux);
		double *auxvec_edge   = FCLAW_ALLOC(double, maux);
		double *flux          = FCLAW_ALLOC(double, meqn);     /* f(qr) - f(ql) = amdq+apdq */

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



	int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
	double* work = FCLAW_ALLOC(double,mwork);

	int size = meqn*(mx+2*mbc)*(my+2*mbc);
	double* fp = FCLAW_ALLOC(double,size);
	double* fm = FCLAW_ALLOC(double,size);
	double* gp = FCLAW_ALLOC(double,size);
	double* gm = FCLAW_ALLOC(double,size);

	int ierror = 0;

	if (claw46_vt->flux2 == NULL)
	{
		claw46_vt->flux2 = (clawpack_options->use_fwaves != 0) ? &CLAWPACK46_FLUX2FW : 
		                       &CLAWPACK46_FLUX2;	
	}

	/* NOTE: qold will be overwritten in this step */
	CLAWPACK46_SET_BLOCK(&blockno);
	CLAWPACK46_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, clawpack_options->method,
						  clawpack_options->mthlim, &clawpack_options->mcapa,
						  &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
						  work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
						  claw46_vt->fort_rpn2, claw46_vt->fort_rpt2,
						  claw46_vt->fort_rpn2fw, claw46_vt->fort_rpt2fw,
						  claw46_vt->flux2,
						  block_corner_count, &ierror, &clawpack_options->use_fwaves);
	CLAWPACK46_UNSET_BLOCK();

	FCLAW_ASSERT(ierror == 0);

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


	FCLAW_FREE(fp);
	FCLAW_FREE(fm);
	FCLAW_FREE(gp);
	FCLAW_FREE(gm);
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
    fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt(glob);

    const fc2d_clawpack46_options_t* clawpack_options;
    clawpack_options = fc2d_clawpack46_get_options(glob);

    if (claw46_vt->b4step2 != NULL)
    {
        fclaw_timer_start_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_B4STEP2]);               
        claw46_vt->b4step2(glob,
                           patch,
                           blockno,
                           patchno,t,dt);

        fclaw_timer_stop_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_B4STEP2]);               
    }

    fclaw_timer_start_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_STEP2]);       

    double maxcfl = clawpack46_step2(glob,
                                     patch,
                                     blockno,
                                     patchno,t,dt);

    fclaw_timer_stop_threadsafe(&glob->timers[FCLAW_TIMER_ADVANCE_STEP2]);       

    if (clawpack_options->src_term > 0 && claw46_vt->src2 != NULL)
    {
        claw46_vt->src2(glob,
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
	const fc2d_clawpack46_options_t* clawpack_options;
	clawpack_options = fc2d_clawpack46_get_options(glob);

	if (clawpack_options->ascii_out != 0)
	{
		fclaw_clawpatch_output_ascii(glob,iframe);
	}

	if (clawpack_options->vtk_out != 0)
	{
		fclaw_clawpatch_output_vtk(glob,iframe);
	}

}



/* ---------------------------------- Virtual table  ---------------------------------- */

static
fc2d_clawpack46_vtable_t* clawpack46_vt_new()
{
    return (fc2d_clawpack46_vtable_t*) FCLAW_ALLOC_ZERO (fc2d_clawpack46_vtable_t, 1);
}

static
void clawpack46_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

void fc2d_clawpack46_solver_initialize(fclaw_global_t* glob)
{
	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
	fc2d_clawpack46_options_t* clawopt = fc2d_clawpack46_get_options(glob);
    clawopt->method[6] = clawpatch_opt->maux;

	if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack46 : bad maux/mcapa combination\n");
        exit(FCLAW_EXIT_ERROR);
    }

	int claw_version = 4;
	fclaw2d_clawpatch_vtable_initialize(glob, claw_version);

	fclaw_vtable_t*                	 fc_vt = fclaw_vt(glob);
	fclaw2d_patch_vtable_t*          patch_vt = fclaw2d_patch_vt(glob);  
    fclaw_clawpatch_vtable_t*      clawpatch_vt = fclaw_clawpatch_vt(glob);

	fc2d_clawpack46_vtable_t*  claw46_vt = clawpack46_vt_new();

	/* ForestClaw vtable items */
	fc_vt->output_frame      = clawpack46_output;
	fc_vt->problem_setup     = clawpack46_setprob;    

	/* These could be over-written by user specific settings */
	patch_vt->initialize                     = clawpack46_qinit;
	patch_vt->setup                          = clawpack46_setaux;  
	patch_vt->physical_bc                    = clawpack46_bc2;
	patch_vt->single_step_update             = clawpack46_update;

	/* Conservation updates (based on Clawpack updates) */
	clawpatch_vt->d2->fort_time_sync_f2c       = CLAWPACK46_FORT_TIME_SYNC_F2C;
	clawpatch_vt->d2->fort_time_sync_samesize  = CLAWPACK46_FORT_TIME_SYNC_SAMESIZE;

	/* Wrappers so that user can change argument list */
	claw46_vt->b4step2                       = clawpack46_b4step2;
	claw46_vt->src2                          = clawpack46_src2;

	/* Required functions  - error if NULL */
	claw46_vt->fort_bc2       = CLAWPACK46_BC2_DEFAULT;
	claw46_vt->fort_qinit     = NULL;
	claw46_vt->fort_rpn2      = NULL;
	claw46_vt->fort_rpt2      = NULL;

	/* Optional functions - call only if non-NULL */
	claw46_vt->fort_setprob   = NULL;
	claw46_vt->fort_setaux    = NULL;
	claw46_vt->fort_b4step2   = NULL;
	claw46_vt->fort_src2      = NULL;

	claw46_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fc2d_clawpack46") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fc2d_clawpack46", claw46_vt, clawpack46_vt_destroy);
}


/* ----------------------------- User access to solver functions --------------------------- */


/* These are here in case the user wants to call Clawpack routines directly */

fc2d_clawpack46_vtable_t* fc2d_clawpack46_vt(fclaw_global_t* glob)
{
	fc2d_clawpack46_vtable_t* claw46_vt = (fc2d_clawpack46_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fc2d_clawpack46");
	FCLAW_ASSERT(claw46_vt != NULL);
	FCLAW_ASSERT(claw46_vt->is_set != 0);
	return claw46_vt;
}

/* This should only be called when a new fclaw_clawpatch_t is created. */
void fc2d_clawpack46_set_capacity(fclaw_global_t *glob,
								  fclaw_patch_t *patch,
								  int blockno,
								  int patchno)
{
	int mx,my,mbc;
	double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	double *area = fclaw2d_clawpatch_get_area(glob,patch);

	int maux;
	double *aux;
	fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

	fc2d_clawpack46_options_t *clawopt = fc2d_clawpack46_get_options(glob);
	int mcapa = clawopt->mcapa;
	FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

	CLAWPACK46_SET_CAPACITY(&mx,&my,&mbc,&dx,&dy,area,&mcapa,
							&maux,aux);
}



/* -------------------------- Public interface to Clawpack wrappers --------------------*/

void fc2d_clawpack46_setprob(fclaw_global_t *glob)
{
	clawpack46_setprob(glob);
}

/* This can be set to claw46_vt->src2 */
void fc2d_clawpack46_src2(fclaw_global_t* glob,
						  fclaw_patch_t *patch,
						  int blockno,
						  int patchno,
						  double t,
						  double dt)
{
	clawpack46_src2(glob,patch,blockno,blockno,t,dt);
}


void fc2d_clawpack46_setaux(fclaw_global_t *glob,
							fclaw_patch_t *patch,
							int blockno,
							int patchno)
{
	clawpack46_setaux(glob,patch,blockno,patchno);
}


void fc2d_clawpack46_qinit(fclaw_global_t *glob,
						   fclaw_patch_t *patch,
						   int blockno,
						   int patchno)
{
	clawpack46_qinit(glob,patch,blockno,patchno);
}

void fc2d_clawpack46_b4step2(fclaw_global_t* glob,
							 fclaw_patch_t *patch,
							 int blockno,
							 int patchno,
							 double t,
							 double dt)
{
	clawpack46_b4step2(glob,patch,blockno,patchno,t,dt);
}

void fc2d_clawpack46_bc2(fclaw_global_t *glob,
						 fclaw_patch_t *patch,
						 int blockno,
						 int patchno,
						 double t,
						 double dt,
						 int intersects_bc[],
						 int time_interp)
{
	clawpack46_bc2(glob,patch,blockno,patchno,t,dt,
				   intersects_bc,time_interp);
}














