/*
  Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "teton_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fc2d_geoclaw.h>
#include <fc2d_geoclaw_options.h>

static
fclaw_domain_t* create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    int mi,mj,a,b;

    mi = fclaw_opt->mi;
    mj = fclaw_opt->mj;
    a = 0; /* non-periodic */
    b = 0;

    /* Rectangular brick domain */
    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
    cont = fclaw2d_map_new_nomap_brick(brick);

    domain = fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    return domain;
}



static
void run_program(fclaw_global_t* glob)
{

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    fc2d_geoclaw_solver_initialize(glob);

    teton_link_solvers(glob);

    fc2d_geoclaw_module_setup(glob);


    /* ---------------------------------------------------------------
       Initialize, run and finalize
       --------------------------------------------------------------- */
    fclaw2d_initialize(glob);
    fclaw2d_run(glob);

    fclaw2d_finalize(glob);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    fclaw_options_t             *fclaw_opt;
    fclaw_clawpatch_options_t *clawpatch_opt;
    fc2d_geoclaw_options_t      *geoclaw_opt;

    fclaw_global_t            *glob;
    fclaw_domain_t            *domain;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    fclaw_opt       =             fclaw_options_register(app,  NULL,       "fclaw_options.ini");
    clawpatch_opt   = fclaw_clawpatch_options_register_2d(app, "clawpatch", "fclaw_options.ini");
    geoclaw_opt     =      fc2d_geoclaw_options_register(app, "geoclaw",   "fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!vexit)
    {
        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw_global_new();
        fclaw_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw2d_options_store           (glob, fclaw_opt);
        fclaw_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_geoclaw_options_store      (glob, geoclaw_opt);

        run_program(glob);
        
        fclaw_global_destroy(glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
