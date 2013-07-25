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

#include "amr_forestclaw.H"
#include "amr_manyclaw.H"
#include "simple_user.H"

#include <manyclaw/manyclaw.h>


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* -----------------------------------------------------------
   Routines below call standard Clawpack routines indirectly through
   calls to ManyClaw interfaces, found in
             forestclaw/src/solvers/amr_manyclaw/amr_manyclaw.cpp

   There are three groups of function pointers that need to be set somewhere.
   fclaw2d_solver_functions_t *sf = get_solver_functions(domain);
   fclaw2d_solver_functions_t *rf = get_regrid_functions(domain);
   fclaw2d_solver_functions_t *of = get_output_functions(domain);

   Solver functions (must all be defined).
   -------------------------------------------------------------------
   sf->f_patch_setup(...)                      --> patch_setup(...) (call setaux, etc)
   sf->f_patch_initialize(...)                 --> patch_initialize(...)  (call qinit, etc)
   sf->f_patch_physical_bc(...)                --> patch_physical_bc(...) (call bc2, etc)
   sf->f_patch_single_step_update(...)         --> patch_single_step_update(...)
                                                   (call b4step2, step2, src2, etc)

   Manyclaw versions of solver functions that can be called
   --------------------------------------------------------
   amr_manyclaw_setprob(...)  --> setprob_()
   amr_manyclaw_setaux(...)   --> setaux_(...)
   amr_manyclaw_qinit(...)    --> qinit_(...)
   amr_manyclaw_bc2(...)      --> bc2_(...)
   amr_manyclaw_b4step2(...)  --> b4step2_(...)
   amr_manyclaw_step2(...)    --> step2_(...) (the amrclaw version;  doesn't update)
   amr_manyclaw_src2(...)     --> src2_(...)

   The routines below are linked to ForestClaw routines with a call to
   'simple_link_solvers'



   Tagging for refinement and coarsening
   -----------------------------------------


   This routine is linked to a problem setup routine using 'link_problem_setup'.
   It is not done above because this is not technically part of a 'solver'. But the
   user could call a solver dependent function here.

   User defined 'simple' routine            linked to ForestClaw as :
   -------------------------------------------------------------------
   simple_problem_setup(...)             --> problem_setup(...) (call setprob, etc)


   --------------------------------------------------------------------------------- */

void simple_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);

    sf->use_single_step_update = fclaw_true;
    sf->use_mol_update = fclaw_false;

    sf->f_patch_setup              = &simple_patch_setup;
    sf->f_patch_initialize         = &simple_patch_initialize;
    sf->f_patch_physical_bc        = &simple_patch_physical_bc;
    sf->f_patch_single_step_update = &simple_patch_single_step_update;

    amr_manyclaw_link_to_clawpatch();
}


void simple_problem_setup(fclaw2d_domain_t* domain)
{
    /* This is called once at the start of the run */
    amr_manyclaw_setprob(domain);
}


void simple_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    /* This is called once when a new patch is created. */
    manyclaw_set_solver(domain,this_patch,this_block_idx,this_patch_idx);
    amr_manyclaw_setaux(domain,this_patch,this_block_idx,this_patch_idx);
    manyclaw_set_riemann_solvers(this_patch,advection_rp_grid_eval_tbb,
                                 updater_first_order_dimensional_splitting);


}



void simple_patch_initialize(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    /* This is called once for each patch in the initial grid layout */
    amr_manyclaw_qinit(domain,this_patch,this_block_idx,this_patch_idx);
}


void simple_patch_physical_bc(fclaw2d_domain *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double dt,
                              fclaw_bool intersects_bc[],
                              fclaw_bool time_interp)
{
    /* This is called everytime a patch needs physical boundary conditions */
    amr_manyclaw_bc2(domain,this_patch,this_block_idx,this_patch_idx,
                     t,dt,intersects_bc,time_interp);
}


double simple_patch_single_step_update(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double dt)
{
    /* This does a single step update on a patch */

    amr_manyclaw_b4step2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);

    double maxcfl = amr_manyclaw_step2(domain,this_patch,this_block_idx,
                                       this_patch_idx,t,dt);
    return maxcfl;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
