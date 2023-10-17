/*
  Copyright (c) 2019-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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
    
#include <fclaw_include_all.h>

#include <fclaw_output.h>
#include <fclaw_diagnostics.h>

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.h>

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_options.h>


static
void create_domain(fclaw_global_t* glob)
{
    /* Mapped, multi-block domain */
    fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
 
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;

    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    fclaw_domain_t *domain = fclaw_domain_new_2d_brick(glob->mpicomm, mi,mj,a,b, fclaw_opt->minlevel);

    /* Map unit square to disk using mapc2m_disk.f */
    fclaw2d_map_context_t *brick = fclaw2d_map_new_brick(domain, mi, mj, a, b);
    fclaw2d_map_context_t *cont = fclaw2d_map_new_nomap_brick(brick);

    fclaw_global_store_domain(glob, domain);
    fclaw2d_map_store(glob, cont);

    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
}

static
void run_program(fclaw_global_t* glob)
{
    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    /* Test thunderegg solver */
    fc2d_thunderegg_solver_initialize(glob);

    /* set up elliptic solver to use the thunderegg solver */
    heat_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */

    /* Set up initial conditions */
    fclaw_initialize(glob);

    heat_run(glob);

    /* ---------------------------------------------------------------
       Finalize
       --------------------------------------------------------------- */
    fclaw_finalize(glob);
}

int
main (int argc, char **argv)
{
    /* Initialize application */
    fclaw_app_t *app = fclaw_app_new (&argc, &argv, NULL);

    /* Options */
    fclaw_options_t             *fclaw_opt;
    fclaw_clawpatch_options_t *clawpatch_opt;
    fc2d_thunderegg_options_t    *mg_opt;
    heat_options_t              *user_opt;

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt =   fclaw_clawpatch_2d_options_register(app, "clawpatch",  "fclaw_options.ini");
    mg_opt =            fc2d_thunderegg_options_register(app, "thunderegg", "fclaw_options.ini");
    user_opt =                     heat_options_register(app,               "fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    int first_arg;
    fclaw_exit_type_t vexit;
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!vexit)
    {
        /* Options have been checked and are valid */
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);
    
        /* Create global structure which stores the domain, timers, etc */
        fclaw_global_t *glob = fclaw_global_new_comm(mpicomm, size, rank);

        /* Store option packages in glob */
        fclaw_options_store           (glob, fclaw_opt);
        fclaw_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_thunderegg_options_store    (glob, mg_opt);
        heat_options_store            (glob, user_opt);

        create_domain(glob);

        run_program(glob);

        fclaw_global_destroy(glob);        
    }
    
    fclaw_app_destroy (app);

    return 0;
}
