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

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"
#include "fc2d_thunderegg_physical_bc.h"
#include "fc2d_thunderegg_fort.h"

#include <fclaw_pointer_map.h>

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>

#include <fclaw2d_patch.h>
#include <fclaw_global.h>
#include <fclaw_vtable.h>
#include <fclaw_output.h>

#include <fclaw_domain.h>

#include "operators/fc2d_thunderegg_fivepoint.h"
#include "operators/fc2d_thunderegg_heat.h"
#include "operators/fc2d_thunderegg_starpatch.h"
#include "operators/fc2d_thunderegg_varpoisson.h"



/* --------------------- ThunderEgg solver (required) ------------------------- */

static
void thunderegg_setup_solver(fclaw_global_t *glob)
{
	//fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);
}


static
void thunderegg_rhs(fclaw_global_t *glob,
                   fclaw_patch_t *patch,
                   int blockno,
                   int patchno)
{
	fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);

	FCLAW_ASSERT(mg_vt->fort_rhs != NULL); /* Must be initialized */

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
	fclaw_clawpatch_rhs_data(glob,patch,&rhs,&mfields);
	FCLAW_ASSERT(mfields == 1);

	/* Compute right hand side */
	mg_vt->fort_rhs(&blockno,&mbc,&mx,&my,&mfields,
                    &xlower,&ylower,&dx,&dy,rhs);
}

static
void thunderegg_solve(fclaw_global_t* glob)
{
    // Apply non-homogeneous boundary conditions 
    fc2d_thunderegg_physical_bc(glob);

    fc2d_thunderegg_vtable_t  *mg_vt  = fc2d_thunderegg_vt(glob);  
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);

    /* Should the operators be part of the thunderegg library? Yes, for now, at least */
    switch (mg_opt->patch_operator)
    {
        case STARPATCH:
            mg_vt->patch_operator = fc2d_thunderegg_starpatch_solve;
            break;
        case FIVEPOINT:
            mg_vt->patch_operator = fc2d_thunderegg_fivepoint_solve;
            break;
        case VARPOISSON:
            mg_vt->patch_operator = fc2d_thunderegg_varpoisson_solve;
            break;
        case HEAT:
            mg_vt->patch_operator = fc2d_thunderegg_heat_solve;
            break;
#if 0
        case USER_OPERATOR:
            if (mg_vt->patch_operator == NULL)
            {
                fclaw_global_essentialf("thunderegg_solve : User specified operator not set\n");
                exit(0);
            }
#endif            
        default:
            break;
            /* user has specified something, hopefully */
    }
    
    FCLAW_ASSERT(mg_vt->patch_operator != NULL);

    mg_vt->patch_operator(glob);
}

/* ---------------------------------- Output functions -------------------------------- */

static
void thunderegg_output(fclaw_global_t *glob, int iframe)
{
	const fc2d_thunderegg_options_t* mg_options;
	mg_options = fc2d_thunderegg_get_options(glob);

	if (mg_options->ascii_out != 0)
	{
		fclaw_clawpatch_output_ascii(glob,iframe);
	}

	if (mg_options->vtk_out != 0)
	{
		fclaw_clawpatch_output_vtk(glob,iframe);
	}
}



/* ------------------------------ Virtual functions  ---------------------------------- */

static
fc2d_thunderegg_vtable_t* thunderegg_vt_new()
{
    return (fc2d_thunderegg_vtable_t*) FCLAW_ALLOC_ZERO (fc2d_thunderegg_vtable_t, 1);
}

static
void thunderegg_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

void fc2d_thunderegg_solver_initialize(fclaw_global_t* glob)
{
	int claw_version = 4; /* solution data is organized as (i,j,m) */
	fclaw2d_clawpatch_vtable_initialize(glob, claw_version);


    //fclaw_clawpatch_vtable_t*      clawpatch_vt = fclaw_clawpatch_vt(glob);

	/* ForestClaw vtable items */
	fclaw_vtable_t*   fc_vt = fclaw_vt(glob);
	fc_vt->output_frame      = thunderegg_output;

	/* These could be over-written by user specific settings */
	fclaw2d_patch_vtable_t*   patch_vt = fclaw2d_patch_vt(glob);  
	patch_vt->rhs            = thunderegg_rhs;  /* Calls FORTRAN routine */
	patch_vt->setup          = NULL;
    
    fclaw_elliptic_vtable_t *elliptic_vt = fclaw_elliptic_vt(glob);
    elliptic_vt->setup = thunderegg_setup_solver;
    elliptic_vt->solve = thunderegg_solve;    
    elliptic_vt->apply_bc = fc2d_thunderegg_physical_bc;

	fc2d_thunderegg_vtable_t*  mg_vt = thunderegg_vt_new();	
    mg_vt->fort_apply_bc = &THUNDEREGG_FORT_APPLY_BC_DEFAULT;
    mg_vt->fort_eval_bc  = &THUNDEREGG_FORT_EVAL_BC_DEFAULT;

#if 0
    /* Operator is specified in solve routine, above*/
    mg_vt->patch_operator = NULL;
#endif    

	mg_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fc2d_thunderegg") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fc2d_thunderegg", mg_vt, thunderegg_vt_destroy);
}


/* ----------------------------- User access to solver functions --------------------------- */

fc2d_thunderegg_vtable_t* fc2d_thunderegg_vt(fclaw_global_t* glob)
{
	fc2d_thunderegg_vtable_t* thunderegg_vt = (fc2d_thunderegg_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fc2d_thunderegg");
	FCLAW_ASSERT(thunderegg_vt != NULL);
	FCLAW_ASSERT(thunderegg_vt->is_set != 0);
	return thunderegg_vt;
}





