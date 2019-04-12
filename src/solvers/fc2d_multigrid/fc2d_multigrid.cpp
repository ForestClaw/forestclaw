/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "fc2d_multigrid.h"
#include "fc2d_multigrid_options.h"

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>

#if 0
#include <fclaw2d_options.h>
#include <fclaw2d_clawpatch_options.h>
#endif


#include <fclaw2d_domain.h>
#include <p4est_bits.h>
#include <p4est_wrap.h>



static fc2d_multigrid_vtable_t s_multigrid_vt;

/* --------------------- Multigrid solver (required) ------------------------- */

static
void multigrid_setup_solver(fclaw2d_global_t *glob)
{
	//fc2d_multigrid_vtable_t*  mg_vt = fc2d_multigrid_vt();
}


static
void multigrid_rhs(fclaw2d_global_t *glob,
                   fclaw2d_patch_t *patch,
                   int blockno,
                   int patchno)
{
	int mx,my,mbc, meqn;
	double dx,dy,xlower,ylower;
	double *q;

	fc2d_multigrid_vtable_t*  mg_vt = fc2d_multigrid_vt();

	FCLAW_ASSERT(mg_vt->fort_rhs != NULL); /* Must be initialized */

	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);
	FCLAW_ASSERT(meqn == 1);

	/* Compute right hand side */
	mg_vt->fort_rhs(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q);
}


static
void multigrid_solve(fclaw2d_global_t* glob)
{
	fc2d_multigrid_solve(glob);
}

/* ---------------------------------- Output functions -------------------------------- */

static
void multigrid_output(fclaw2d_global_t *glob, int iframe)
{
	const fc2d_multigrid_options_t* mg_options;
	mg_options = fc2d_multigrid_get_options(glob);

	if (mg_options->ascii_out != 0)
	{
		fclaw2d_clawpatch_output_ascii(glob,iframe);
	}

	if (mg_options->vtk_out != 0)
	{
		fclaw2d_clawpatch_output_vtk(glob,iframe);
	}
}



/* ------------------------------ Virtual functions  ---------------------------------- */

static
fc2d_multigrid_vtable_t* multigrid_vt_init()
{
	FCLAW_ASSERT(s_multigrid_vt.is_set == 0);
	return &s_multigrid_vt;
}

void fc2d_multigrid_solver_initialize()
{
	int claw_version = 4; /* solution data is organized as (i,j,m) */
	fclaw2d_clawpatch_vtable_initialize(claw_version);


	fc2d_multigrid_vtable_t*  mg_vt = multigrid_vt_init();	

    //fclaw2d_clawpatch_vtable_t*      clawpatch_vt = fclaw2d_clawpatch_vt();


	/* ForestClaw vtable items */
	fclaw2d_vtable_t*   fclaw_vt = fclaw2d_vt();
	fclaw_vt->output_frame      = multigrid_output;

	/* These could be over-written by user specific settings */
	fclaw2d_patch_vtable_t*   patch_vt = fclaw2d_patch_vt();  
	patch_vt->rhs            = multigrid_rhs;  /* Calls FORTRAN routine */
	patch_vt->setup          = NULL;

    
    fclaw2d_elliptic_vtable_t *elliptic_vt = fclaw2d_elliptic_vt();
    elliptic_vt->setup = multigrid_setup_solver;
    elliptic_vt->solve = multigrid_solve;
    // RHS is set using default RHS in elliptic_solve.

	mg_vt->is_set = 1;
}


/* ----------------------------- User access to solver functions --------------------------- */

fc2d_multigrid_vtable_t* fc2d_multigrid_vt()
{
	FCLAW_ASSERT(s_multigrid_vt.is_set != 0);
	return &s_multigrid_vt;
}





