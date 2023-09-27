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

#include "heat_user.h"
    
#include <fclaw_filesystem.h>

#include <fclaw2d_include_all.h>

#include <fclaw2d_output.h>
#include <fclaw2d_global.h>
#include <fclaw2d_diagnostics.h>

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch.h>

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_options.h>


void heat_create_domain(fclaw2d_global_t* glob)
{
    /* Mapped, multi-block domain */
    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
 
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;

    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    fclaw2d_domain_t *domain = fclaw2d_domain_new_brick(glob->mpicomm, mi,mj,a,b, fclaw_opt->minlevel);

    /* Map unit square to disk using mapc2m_disk.f */
    fclaw2d_map_context_t *brick = fclaw2d_map_new_brick(domain, mi, mj, a, b);
    fclaw2d_map_context_t *cont = fclaw2d_map_new_nomap_brick(brick);

    fclaw2d_global_store_domain(glob, domain);
    fclaw2d_global_store_map(glob, cont);

    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
}

void heat_run_program(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);


    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Test thunderegg solver */
    fc2d_thunderegg_solver_initialize(glob);

    /* set up elliptic solver to use the thunderegg solver */
    heat_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */

    /* Set up initial conditions */
    fclaw2d_initialize(glob);

    heat_run(glob);

    /* ---------------------------------------------------------------
       Finalize
       --------------------------------------------------------------- */
    fclaw2d_finalize(glob);

    fclaw2d_clear_global_context(glob);
}