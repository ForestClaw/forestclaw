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

#include <p4est_base.h>
#include "amr_utils.H"
#include "amr_single_step.H"
#include "amr_mol.H"
#include "amr_solver_typedefs.H"
#include "fclaw2d_solvers.H"

void amr_dummy_patch_setup(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    /* This is called if there is nothing to do to set up a new patch */
}


void initialize_solver_functions(fclaw2d_solver_functions_t* solver_functions)
{
    fclaw2d_solver_functions_t *sf = solver_functions;

    sf->f_patch_setup              = NULL;
    sf->f_patch_initialize         = NULL;
    sf->f_patch_physical_bc        = NULL;
    sf->f_level_single_step        = &amr_level_single_step_update;  /* default? */
    sf->f_patch_single_step_update = NULL;
    sf->f_level_ode_solver         = NULL;
    sf->f_patch_ode_solver_rhs     = NULL;
}


void copy_solver_functions(fclaw2d_solver_functions_t* old_solver_functions,
                           fclaw2d_solver_functions_t* new_solver_functions)
{

    fclaw2d_solver_functions_t *oldsf = old_solver_functions;
    fclaw2d_solver_functions_t *newsf = new_solver_functions;

    newsf->f_patch_setup              = oldsf->f_patch_setup;
    newsf->f_patch_initialize         = oldsf->f_patch_initialize;
    newsf->f_patch_physical_bc        = oldsf->f_patch_physical_bc;
    newsf->f_level_single_step        = oldsf->f_level_single_step;
    newsf->f_patch_single_step_update = oldsf->f_patch_single_step_update;
    newsf->f_level_ode_solver         = oldsf->f_level_ode_solver;
    newsf->f_patch_ode_solver_rhs     = oldsf->f_patch_ode_solver_rhs;
}
