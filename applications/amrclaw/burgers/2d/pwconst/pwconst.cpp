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

#include <amr_single_step.h>
#include <fc2d_clawpack46.H>
#include <fclaw2d_map.h>
#include <p4est_connectivity.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw2d_map_query.h>

#include "pwconst_user.H"

void run_program(fclaw_app_t* app, amr_options_t* gparms,
                 fc2d_clawpack46_options_t* clawpack_options)
{
    sc_MPI_Comm            mpicomm;

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t	     *domain;
    fclaw2d_map_context_t    *cont = NULL;

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

    /* Use [ax,bx]x[ay,by] */
    conn = p4est_connectivity_new_unitsquare();
    cont = fclaw2d_map_new_nomap();

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

    /* ---------------------------------------------------------- */
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_INFO);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    init_domain_data(domain);

    set_domain_parms(domain,gparms);
    fc2d_clawpack46_set_options(domain,clawpack_options);

    pwconst_link_solvers(domain);

    link_regrid_functions(domain,pwconst_patch_tag4refinement,
                          pwconst_patch_tag4coarsening);

    amrinit(&domain);
    amrrun(&domain);
    amrreset(&domain);

    fclaw2d_map_destroy(cont);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                  *options;
    amr_options_t                 samr_options, *gparms = &samr_options;
    fc2d_clawpack46_options_t     sclawpack_options, *clawpack_options = &sclawpack_options;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);
    options = fclaw_app_get_options (app);

    fclaw_options_register_general (app, "fclaw_options.ini", gparms);
    fc2d_clawpack46_options_register(app,NULL,clawpack_options);

    /* Read configuration file(s) */
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* -------------------------------------------------------------
       - Run program
       ------------------------------------------------------------- */
    if (!retval & !vexit)
    {
        run_program(app, gparms, clawpack_options);
    }

    fclaw_app_destroy (app);

    return 0;
}
