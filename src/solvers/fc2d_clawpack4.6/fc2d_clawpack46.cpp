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

#include "fc2d_clawpack46.h"
#include "fc2d_clawpack46_options.h"
#include "fc2d_clawpack46_fort.h"

#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_clawpatch.h>

#include <fclaw2d_clawpatch_diagnostics.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fclaw2d_clawpatch_conservation.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_options.h>
#include <fclaw2d_defs.h>


static fc2d_clawpack46_vtable_t s_clawpack46_vt;

/* --------------------- Clawpack solver functions (required) ------------------------- */

static
void clawpack46_setprob(fclaw2d_global_t *glob)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();
	if (claw46_vt->fort_setprob != NULL)
	{
		claw46_vt->fort_setprob();
	}
}


static
void clawpack46_qinit(fclaw2d_global_t *glob,
					  fclaw2d_patch_t *this_patch,
					  int this_block_idx,
					  int this_patch_idx)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();

	FCLAW_ASSERT(claw46_vt->fort_qinit != NULL); /* Must be initialized */
	int mx,my,mbc,meqn,maux,maxmx,maxmy;
	double dx,dy,xlower,ylower;
	double *q, *aux;

	fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

	maxmx = mx;
	maxmy = my;

	/* Call to classic Clawpack 'qinit' routine.  This must be user defined */
	CLAWPACK46_SET_BLOCK(&this_block_idx);
	claw46_vt->fort_qinit(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
						  &maux,aux);
	CLAWPACK46_UNSET_BLOCK();
}


static
void clawpack46_bc2(fclaw2d_global_t *glob,
					fclaw2d_patch_t *this_patch,
					int this_block_idx,
					int this_patch_idx,
					double t,
					double dt,
					int intersects_phys_bdry[],
					int time_interp)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();

	fc2d_clawpack46_options_t *clawpack_options = fc2d_clawpack46_get_options(glob);

	FCLAW_ASSERT(claw46_vt->fort_bc2 != NULL);

	int mx,my,mbc,meqn,maux,maxmx,maxmy;
	double xlower,ylower,dx,dy;
	double *aux,*q;

	fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

	maxmx = mx;
	maxmy = my;

	int *block_mthbc = clawpack_options->mthbc;

	/* Set a local copy of mthbc that can be used for a patch. */
	int mthbc[4];
	for(int i = 0; i < 4; i++)
	{
		if (intersects_phys_bdry[i])
		{
			mthbc[i] = block_mthbc[i];
		}
		else
		{
			mthbc[i] = -1;
		}
	}

	/*
	  We may be imposing boundary conditions on time-interpolated data;
	  and is being done just to help with fine grid interpolation.
	  In this case, this boundary condition won't be used to update
	  anything
	*/
	fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&q,&meqn);

	CLAWPACK46_SET_BLOCK(&this_block_idx);
	claw46_vt->fort_bc2(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,
						&dx,&dy,q,&maux,aux,&t,&dt,mthbc);
	CLAWPACK46_UNSET_BLOCK();
}


static
void clawpack46_b4step2(fclaw2d_global_t *glob,
						fclaw2d_patch_t *this_patch,
						int this_block_idx,
						int this_patch_idx,
						double t, double dt)

{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();

	int mx,my,mbc,meqn, maux,maxmx,maxmy;
	double xlower,ylower,dx,dy;
	double *aux,*q;

	if (claw46_vt->fort_b4step2 == NULL)
	{
		return;
	}

	fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

	maxmx = mx;
	maxmy = my;

	CLAWPACK46_SET_BLOCK(&this_block_idx);
	claw46_vt->fort_b4step2(&maxmx,&maxmy,&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
							&dx,&dy,&t,&dt,&maux,aux);
	CLAWPACK46_UNSET_BLOCK();
}

static
void clawpack46_src2(fclaw2d_global_t *glob,
					 fclaw2d_patch_t *this_patch,
					 int this_block_idx,
					 int this_patch_idx,
					 double t,
					 double dt)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();

	int mx,my,mbc,meqn, maux,maxmx,maxmy;
	double xlower,ylower,dx,dy;
	double *aux,*q;

	if (claw46_vt->fort_src2 == NULL)
	{
		/* User has not set a fortran routine */
		return;
	}

	fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

	maxmx = mx;
	maxmy = my;

	CLAWPACK46_SET_BLOCK(&this_block_idx);
	claw46_vt->fort_src2(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,
						 &dx,&dy,q,&maux,aux,&t,&dt);
	CLAWPACK46_UNSET_BLOCK();
}


/* This can be used as a value for patch_vt->patch_setup */
static
void clawpack46_setaux(fclaw2d_global_t *glob,
					   fclaw2d_patch_t *this_patch,
					   int this_block_idx,
					   int this_patch_idx)
	{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();

	if (claw46_vt->fort_setaux == NULL)
	{
		return;
	}

	if (fclaw2d_patch_is_ghost(this_patch))
	{
		/* This is going to be removed at some point */
		return;
	}

	int mx,my,mbc,maux,maxmx,maxmy;
	double xlower,ylower,dx,dy;
	double *aux;

	fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

	maxmx = mx;
	maxmy = my;

	CLAWPACK46_SET_BLOCK(&this_block_idx);
	claw46_vt->fort_setaux(&maxmx,&maxmy,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
						   &maux,aux);
	CLAWPACK46_UNSET_BLOCK();
}

/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
static
double clawpack46_step2(fclaw2d_global_t *glob,
						fclaw2d_patch_t *this_patch,
						int this_block_idx,
						int this_patch_idx,
						double t,
						double dt)
{
	fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();
	const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
	const fc2d_clawpack46_options_t* clawpack_options;
	int level;
	double *qold, *aux;
	int mx, my, meqn, maux, mbc;
	double xlower, ylower, dx,dy;

	clawpack_options = fc2d_clawpack46_get_options(glob);

	FCLAW_ASSERT(claw46_vt->fort_rpn2 != NULL);

	if (clawpack_options->order[1] > 0)
	{
		FCLAW_ASSERT(claw46_vt->fort_rpt2 != NULL);
	}


	level = this_patch->level;

	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

	fclaw2d_clawpatch_save_current_step(glob, this_patch);

	fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);


	int mwaves = clawpack_options->mwaves;

	int maxm = fmax(mx,my);

	double cflgrid = 0.0;

	fclaw2d_clawpatch_registers_t* cr = 
		  fclaw2d_clawpatch_get_registers(glob,this_patch);

	/* Evaluate fluxes needed in correction terms */
	if (fclaw_opt->time_sync && fclaw_opt->flux_correction)
	{
		FCLAW_ASSERT(claw46_vt->fort_rpn2_cons != NULL);
		double *qvec   = FCLAW_ALLOC(double, meqn);
		double *auxvec_center = FCLAW_ALLOC(double, maux);
		double *auxvec_edge = FCLAW_ALLOC(double, maux);
		double *flux   = FCLAW_ALLOC(double, meqn);     /* f(qr) - f(ql) = amdq+apdq */

		CLAWPACK46_TIME_SYNC_STORE_FLUX(&mx,&my,&mbc,&meqn,&maux,
		                                &this_block_idx,&this_patch_idx, &dt,
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
	double* work = new double[mwork];

	int size = meqn*(mx+2*mbc)*(my+2*mbc);
	double* fp = new double[size];
	double* fm = new double[size];
	double* gp = new double[size];
	double* gm = new double[size];

	int ierror = 0;
	int* block_corner_count = fclaw2d_patch_block_corner_count(glob,this_patch);
	clawpack46_fort_flux2_t flux2 = clawpack_options->use_fwaves ?
									 CLAWPACK46_FLUX2FW : CLAWPACK46_FLUX2;

	/* NOTE: qold will be overwritten in this step */
	CLAWPACK46_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, clawpack_options->method,
						  clawpack_options->mthlim, &clawpack_options->mcapa,
						  &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
						  work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
						  claw46_vt->fort_rpn2, claw46_vt->fort_rpt2,flux2,
						  block_corner_count, &ierror);

	FCLAW_ASSERT(ierror == 0);

	if (fclaw_opt->time_sync && fclaw_opt->fluctuation_correction)
	{

		CLAWPACK46_TIME_SYNC_ACCUMULATE_WAVES(&mx,&my,&mbc,&meqn, &dt, &dx, 
		                                      &dy, &this_patch_idx,
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


	delete [] fp;
	delete [] fm;
	delete [] gp;
	delete [] gm;

	delete [] work;

	return cflgrid;
}

static
double clawpack46_update(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt, 
                         void* user)
{
    fc2d_clawpack46_vtable_t*  claw46_vt = fc2d_clawpack46_vt();

    const fc2d_clawpack46_options_t* clawpack_options;
    clawpack_options = fc2d_clawpack46_get_options(glob);

    if (claw46_vt->b4step2 != NULL)
    {
        fclaw2d_timer_start_threadsafe(&glob->timers[FCLAW2D_TIMER_ADVANCE_B4STEP2]);               
        claw46_vt->b4step2(glob,
                           this_patch,
                           this_block_idx,
                           this_patch_idx,t,dt);

        fclaw2d_timer_stop_threadsafe(&glob->timers[FCLAW2D_TIMER_ADVANCE_B4STEP2]);               
    }

    fclaw2d_timer_start_threadsafe(&glob->timers[FCLAW2D_TIMER_ADVANCE_STEP2]);       

    double maxcfl = clawpack46_step2(glob,
                                     this_patch,
                                     this_block_idx,
                                     this_patch_idx,t,dt);

    fclaw2d_timer_stop_threadsafe(&glob->timers[FCLAW2D_TIMER_ADVANCE_STEP2]);       

    if (clawpack_options->src_term > 0 && claw46_vt->src2 != NULL)
    {
        claw46_vt->src2(glob,
                        this_patch,
                        this_block_idx,
                        this_patch_idx,t,dt);
    }
    return maxcfl;
}


/* ---------------------------------- Output functions -------------------------------- */

static
void clawpack46_output(fclaw2d_global_t *glob, int iframe)
{
	const fc2d_clawpack46_options_t* clawpack_options;
	clawpack_options = fc2d_clawpack46_get_options(glob);

	if (clawpack_options->ascii_out != 0)
	{
		fclaw2d_clawpatch_output_ascii(glob,iframe);
	}

	if (clawpack_options->vtk_out != 0)
	{
		fclaw2d_clawpatch_output_vtk(glob,iframe);
	}

}



/* ---------------------------------- Virtual table  ---------------------------------- */

static
fc2d_clawpack46_vtable_t* clawpack46_vt_init()
{
	FCLAW_ASSERT(s_clawpack46_vt.is_set == 0);
	return &s_clawpack46_vt;
}

void fc2d_clawpack46_solver_initialize()
{
	int claw_version = 4;
	fclaw2d_clawpatch_vtable_initialize(claw_version);

	fclaw2d_vtable_t*                fclaw_vt = fclaw2d_vt();
	fclaw2d_patch_vtable_t*          patch_vt = fclaw2d_patch_vt();  
    fclaw2d_clawpatch_vtable_t*      clawpatch_vt = fclaw2d_clawpatch_vt();

	fc2d_clawpack46_vtable_t*  claw46_vt = clawpack46_vt_init();

	/* ForestClaw vtable items */
	fclaw_vt->output_frame      = clawpack46_output;
	fclaw_vt->problem_setup     = clawpack46_setprob;    

	/* These could be over-written by user specific settings */
	patch_vt->initialize                     = clawpack46_qinit;
	patch_vt->setup                          = clawpack46_setaux;  
	patch_vt->physical_bc                    = clawpack46_bc2;
	patch_vt->single_step_update             = clawpack46_update;

	/* Conservation updates (based on Clawpack updates) */
	clawpatch_vt->fort_time_sync_f2c         = CLAWPACK46_FORT_TIME_SYNC_F2C;
	clawpatch_vt->fort_time_sync_samesize    = CLAWPACK46_FORT_TIME_SYNC_SAMESIZE;

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
}


/* ----------------------------- User access to solver functions --------------------------- */


/* These are here in case the user wants to call Clawpack routines directly */

fc2d_clawpack46_vtable_t* fc2d_clawpack46_vt()
{
	FCLAW_ASSERT(s_clawpack46_vt.is_set != 0);
	return &s_clawpack46_vt;
}

/* This should only be called when a new fclaw2d_clawpatch_t is created. */
void fc2d_clawpack46_set_capacity(fclaw2d_global_t *glob,
								  fclaw2d_patch_t *this_patch,
								  int this_block_idx,
								  int this_patch_idx)
{
	int mx,my,mbc,maux,mcapa;
	double dx,dy,xlower,ylower;
	double *aux, *area;
	fc2d_clawpack46_options_t *clawopt;

	clawopt = fc2d_clawpack46_get_options(glob);
	mcapa = clawopt->mcapa;

	fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	area = fclaw2d_clawpatch_get_area(glob,this_patch);

	fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
	FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

	CLAWPACK46_SET_CAPACITY(&mx,&my,&mbc,&dx,&dy,area,&mcapa,
							&maux,aux);
}



/* -------------------------- Public interface to Clawpack wrappers --------------------*/

/* These are overkill;  it isn't obvious why the user would want these */

void fc2d_clawpack46_setprob(fclaw2d_global_t *glob)
{
	clawpack46_setprob(glob);
}

/* This can be set to claw46_vt->src2 */
void fc2d_clawpack46_src2(fclaw2d_global_t* glob,
						  fclaw2d_patch_t *this_patch,
						  int this_block_idx,
						  int this_patch_idx,
						  double t,
						  double dt)
{
	clawpack46_src2(glob,this_patch,this_block_idx,this_block_idx,t,dt);
}


void fc2d_clawpack46_setaux(fclaw2d_global_t *glob,
							fclaw2d_patch_t *this_patch,
							int this_block_idx,
							int this_patch_idx)
{
	clawpack46_setaux(glob,this_patch,this_block_idx,this_patch_idx);
}


void fc2d_clawpack46_qinit(fclaw2d_global_t *glob,
						   fclaw2d_patch_t *this_patch,
						   int this_block_idx,
						   int this_patch_idx)
{
	clawpack46_qinit(glob,this_patch,this_block_idx,this_patch_idx);
}

void fc2d_clawpack46_b4step2(fclaw2d_global_t* glob,
							 fclaw2d_patch_t *this_patch,
							 int this_block_idx,
							 int this_patch_idx,
							 double t,
							 double dt)
{
	clawpack46_b4step2(glob,this_patch,this_block_idx,this_patch_idx,t,dt);
}

void fc2d_clawpack46_bc2(fclaw2d_global_t *glob,
						 fclaw2d_patch_t *this_patch,
						 int this_block_idx,
						 int this_patch_idx,
						 double t,
						 double dt,
						 int intersects_bc[],
						 int time_interp)
{
	clawpack46_bc2(glob,this_patch,this_block_idx,this_block_idx,t,dt,
				   intersects_bc,time_interp);
}














