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
#include "fclaw2d_solvers.H"
#include "fclaw2d_defs.h"

void amr_dummy_patch_setup(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    /* Do nothing */
}

void amr_dummy_patch_initialize(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{
    /* Do nothing */
}


void amr_dummy_patch_physical_bc(fclaw2d_domain *domain,
                                 fclaw2d_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx,
                                 double t,
                                 double dt,
                                 fclaw_bool *intersects_bc)
{
    /* Do nothing */
}

void amr_dummy_level_ode_solver(int neqn, double q[], double t, double dt)
{
    /* Do nothing */
}

double amr_dummy_patch_ode_solver_rhs(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double *rhs)
{
    /* Do nothing */
    return 0;
}

double amr_dummy_level_single_step(fclaw2d_domain_t *domain,
                                   int level,
                                   double t, double dt)
{
    /* Do nothing */
    return 0;
}

double amr_dummy_patch_single_step_update(fclaw2d_domain_t *domain,
                                          fclaw2d_patch_t *this_patch,
                                          int this_block_idx,
                                          int this_patch_idx,
                                          double t,
                                          double dt)
{
    /* Do nothing */
    return 0;
}


void initialize_solver_functions(fclaw2d_solver_functions_t* solver_functions)
{
    fclaw2d_solver_functions_t *sf = solver_functions;

    sf->use_single_step_update = fclaw_true;
    sf->use_mol_update = fclaw_false;

    sf->f_patch_setup              = &amr_dummy_patch_setup;
    sf->f_patch_initialize         = &amr_dummy_patch_initialize;
    sf->f_patch_physical_bc        = &amr_dummy_patch_physical_bc;

    sf->f_level_single_step        = &amr_level_single_step_update;  /* default? */

    sf->f_patch_single_step_update = &amr_dummy_patch_single_step_update;

    sf->f_level_ode_solver         = &amr_dummy_level_ode_solver;
    sf->f_patch_ode_solver_rhs     = &amr_dummy_patch_ode_solver_rhs;
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
