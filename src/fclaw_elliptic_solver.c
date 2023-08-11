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

#include <fclaw_pointer_map.h>

#include <fclaw_elliptic_solver.h>
#include <fclaw_global.h>

#include <fclaw2d_patch.h>


/* ---------------------------- Setup solver functions -------------------------------- */

static
void elliptic_setup_solver(fclaw_global_t *glob)
{
    fclaw_elliptic_vtable_t* elliptic_vt = fclaw_elliptic_vt(glob);

    /* May not be anything to set up */
    if (elliptic_vt->setup != NULL)
        elliptic_vt->setup(glob);
}

/* -------------------------------- RHS functions -------------------------------- */

static
void cb_elliptic_rhs(fclaw_domain_t *domain,
                     fclaw_patch_t *patch,
                     int blockno,
                     int patchno,
                     void* user)
{
    fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;

    fclaw_patch_vtable_t* patch_vt = fclaw2d_patch_vt(g->glob);

    /* Check that the user has set the right hand side function */
    FCLAW_ASSERT(patch_vt->rhs != NULL);

    patch_vt->rhs(g->glob,patch,blockno, patchno);    
}

static
void elliptic_rhs_default(fclaw_global_t *glob)
{
    /* Set up right hand side by iterating over patches. */
    fclaw_global_iterate_patches (glob, cb_elliptic_rhs, NULL);
}

static
void elliptic_rhs(fclaw_global_t *glob)
{
    fclaw_elliptic_vtable_t* elliptic_vt = fclaw_elliptic_vt(glob);

    FCLAW_ASSERT(elliptic_vt->rhs != NULL);

    /* This is virtualized in case users have some other way of setting up a right
       hand side (using data from a time dependent solution, for example) */

    elliptic_vt->rhs(glob);
}

/* ----------------------------------- Solve functions -------------------------------- */

static
void elliptic_solve(fclaw_global_t *glob)
{
    fclaw_elliptic_vtable_t* elliptic_vt = fclaw_elliptic_vt(glob);

    /* Check that the user has set an elliptic solver */
    FCLAW_ASSERT(elliptic_vt->solve != NULL);
    
    elliptic_vt->solve(glob);
}


/*----------------------------------- Public interface -------------------------------- */

void fclaw_elliptic_solve(fclaw_global_t *glob)
{
    fclaw_domain_t* domain = glob->domain;

    /* Any solver related setup that should happen before each solve */
    elliptic_setup_solver(glob);

    /* Set up right hand side */
    elliptic_rhs(glob);

    /* Pass p4est pointer to elliptic solver, solve and transfer solution back. */
    glob->count_elliptic_grids +=  domain->local_num_patches;
    fclaw_timer_start (&glob->timers[FCLAW_TIMER_ELLIPTIC_SOLVE]);    

    elliptic_solve(glob);
    
    fclaw_timer_stop (&glob->timers[FCLAW_TIMER_ELLIPTIC_SOLVE]);    

}

/*---------------------------- Virtual table functions -------------------------------- */

static
fclaw_elliptic_vtable_t* elliptic_vt_new()
{
    return (fclaw_elliptic_vtable_t*) FCLAW_ALLOC_ZERO (fclaw_elliptic_vtable_t, 1);
}

static
void elliptic_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}


void fclaw_elliptic_vtable_initialize(fclaw_global_t* glob)
{
    fclaw_elliptic_vtable_t *elliptic_vt = elliptic_vt_new();

    elliptic_vt->rhs = elliptic_rhs_default;

    /* User should set the setup, rhs and solve routines */
    elliptic_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fclaw2d_elliptic") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fclaw2d_elliptic", elliptic_vt, elliptic_vt_destroy);
}

/*----------------------------------- Access functions -------------------------------- */

fclaw_elliptic_vtable_t* fclaw_elliptic_vt(fclaw_global_t* glob)
{
	fclaw_elliptic_vtable_t* elliptic_vt = (fclaw_elliptic_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fclaw2d_elliptic");
	FCLAW_ASSERT(elliptic_vt != NULL);
	FCLAW_ASSERT(elliptic_vt->is_set != 0);
	return elliptic_vt;
}


