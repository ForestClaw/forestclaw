/*
Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "fc3d_thunderegg.h"
#include "fc3d_thunderegg_options.h"
#include "fc3d_thunderegg_physical_bc.h"
#include "fc3d_thunderegg_fort.h"

#include <fclaw_pointer_map.h>

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>
#ifdef FCLAW_WITH_HDF5
#include <fclaw_clawpatch_output_hdf5.h>
#endif

#include <fclaw_patch.h>
#include <fclaw_global.h>
#include <fclaw_vtable.h>
#include <fclaw_output.h>

#include <fclaw_domain.h>

//#include "operators/fc3d_thunderegg_fivepoint.h"
//#include "operators/fc3d_thunderegg_heat.h"
//#include "operators/fc3d_thunderegg_starpatch.h"
//#include "operators/fc3d_thunderegg_varpoisson.h"



/* --------------------- ThunderEgg solver (required) ------------------------- */

static
void thunderegg_setup_solver(fclaw_global_t *glob)
{
	//fc3d_thunderegg_vtable_t*  mg_vt = fc3d_thunderegg_vt(glob);
}


static
void thunderegg_rhs(fclaw_global_t *glob,
                   fclaw_patch_t *patch,
                   int blockno,
                   int patchno)
{
	fc3d_thunderegg_vtable_t*  mg_vt = fc3d_thunderegg_vt(glob);

	FCLAW_ASSERT(mg_vt->fort_rhs != NULL); /* Must be initialized */

    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
	fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz,&mbc,
								&xlower,&ylower,&zlower,&dx,&dy,&dz);

    int mfields;
    double *rhs;
	fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);
	FCLAW_ASSERT(mfields == 1);

	/* Compute right hand side */
	mg_vt->fort_rhs(&blockno,&mbc,&mx,&my,&mz,&mfields,
                    &xlower,&ylower,&zlower,&dx,&dy,&dz,rhs);
}

static
void thunderegg_solve(fclaw_global_t* glob)
{
    // Apply non-homogeneous boundary conditions 
    fc3d_thunderegg_physical_bc(glob);

    fc3d_thunderegg_vtable_t  *mg_vt  = fc3d_thunderegg_vt(glob);  

    FCLAW_ASSERT(mg_vt->patch_operator != NULL);

    mg_vt->patch_operator(glob);
}

/* ---------------------------------- Output functions -------------------------------- */

static
void thunderegg_output(fclaw_global_t *glob, int iframe)
{
	const fc3d_thunderegg_options_t* mg_options;
	mg_options = fc3d_thunderegg_get_options(glob);

	if (mg_options->ascii_out != 0)
	{
		fclaw_clawpatch_output_ascii(glob,iframe);
	}

	if (mg_options->vtk_out != 0)
	{
		fclaw_clawpatch_output_vtk(glob,iframe);
	}

#ifdef FCLAW_WITH_HDF5
	if (mg_options->hdf5_out != 0)
	{
		fclaw_clawpatch_output_hdf5(glob,iframe);
	}
#endif
}



/* ------------------------------ Virtual functions  ---------------------------------- */

static
fc3d_thunderegg_vtable_t* thunderegg_vt_new()
{
    return (fc3d_thunderegg_vtable_t*) FCLAW_ALLOC_ZERO (fc3d_thunderegg_vtable_t, 1);
}

static
void thunderegg_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

void fc3d_thunderegg_solver_initialize(fclaw_global_t* glob)
{
	int claw_version = 4; /* solution data is organized as (i,j,m) */
	fclaw_clawpatch_vtable_initialize(glob, claw_version);


    //fclaw3d_clawpatch_vtable_t*      clawpatch_vt = fclaw3d_clawpatch_vt(glob);

	/* ForestClaw vtable items */
	fclaw_vtable_t*   fc_vt = fclaw_vt(glob);
	fc_vt->output_frame      = thunderegg_output;

	/* These could be over-written by user specific settings */
	fclaw_patch_vtable_t*   patch_vt = fclaw_patch_vt(glob);  
	patch_vt->rhs            = thunderegg_rhs;  /* Calls FORTRAN routine */
	patch_vt->setup          = NULL;
    
    fclaw_elliptic_vtable_t *elliptic_vt = fclaw_elliptic_vt(glob);
    elliptic_vt->setup = thunderegg_setup_solver;
    elliptic_vt->solve = thunderegg_solve;    
    elliptic_vt->apply_bc = fc3d_thunderegg_physical_bc;

	fc3d_thunderegg_vtable_t*  mg_vt = thunderegg_vt_new();	
    mg_vt->fort_apply_bc = &FC3D_THUNDEREGG_FORT_APPLY_BC_DEFAULT;
    mg_vt->fort_eval_bc  = &FC3D_THUNDEREGG_FORT_EVAL_BC_DEFAULT;

    mg_vt->patch_operator = NULL;

	mg_vt->is_set = 1;

	fclaw_global_vtable_store(glob, "fc3d_thunderegg", mg_vt, thunderegg_vt_destroy);
}


/* ----------------------------- User access to solver functions --------------------------- */

fc3d_thunderegg_vtable_t* fc3d_thunderegg_vt(fclaw_global_t* glob)
{
	fc3d_thunderegg_vtable_t* thunderegg_vt = (fc3d_thunderegg_vtable_t*) 
		fclaw_global_get_vtable(glob, "fc3d_thunderegg");
	FCLAW_ASSERT(thunderegg_vt != NULL);
	FCLAW_ASSERT(thunderegg_vt->is_set != 0);
	return thunderegg_vt;
}





