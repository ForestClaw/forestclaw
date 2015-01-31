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

#include <fclaw_options.h>
#include <amr_forestclaw.H>
#include <amr_utils.H>

#include <fclaw2d_clawpack.H>
#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>

#include <fclaw_options.h>

#include <p4est_connectivity.h>

#include "swirl_user.H"

void run_program(fclaw_app_t* app, amr_options_t* gparms,
                 clawpack46_options_t* clawpack_options)
{
    sc_MPI_Comm            mpicomm;

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t	     *domain;
    fclaw2d_map_context_t    *cont = NULL;

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

    /* ---------------------------------------------------------------
       Domain geometry
       -------------------------------------------------------------- */

    /* Map unit square to disk using mapc2m_disk.f */
    gparms->manifold = 0;
    conn = p4est_connectivity_new_unitsquare();
    cont = fclaw2d_map_new_nomap();

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    init_domain_data(domain);

    set_domain_parms(domain,gparms);
    set_clawpack46_options(domain,clawpack_options);

    /* ---------------------------------------------------------------
       Define the solver and link in other problem/user specific
       routines
       --------------------------------------------------------------- */

    link_problem_setup(domain,swirl_problem_setup);

    swirl_link_solvers(domain);

    link_regrid_functions(domain,swirl_patch_tag4refinement,
                          swirl_patch_tag4coarsening);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    amrinit(&domain);
    amrrun(&domain);
    amrreset(&domain);

    /* This has to be in this scope */
    fclaw2d_map_destroy(cont);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t             *options;
    amr_options_t            samr_options,      *gparms = &samr_options;
    clawpack46_options_t     sclawpack_options, *clawpack_options = &sclawpack_options;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);
    options = fclaw_app_get_options (app);

    /*  Register options for each package */
    fclaw_app_options_register_general (app, "fclaw_options.ini", gparms);
    clawpack46_app_options_register (app, "fclaw_options.ini", clawpack_options);

    /* Read configuration file(s) and command line, and process options */
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!retval & !vexit)
    {
        run_program(app, gparms, clawpack_options);
    }

    fclaw_app_destroy (app);

    return 0;
}
