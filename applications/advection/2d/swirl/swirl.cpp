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

#include <p4est_connectivity.h>


#include "swirl_user.H"

int
main (int argc, char **argv)
{
    int		          lp;
    sc_MPI_Comm           mpicomm;
    sc_options_t          *options;
    p4est_connectivity_t  *conn = NULL;
    fclaw2d_map_context_t *cont = NULL;
    fclaw2d_domain_t	  *domain;

    amr_options_t  samr_options, *gparms = &samr_options;
    fclaw2d_clawpack_parms_t  sclawparms, *clawpack_parms = &sclawparms;

    int retval;

    lp = SC_LP_PRODUCTION;
    mpicomm = sc_MPI_COMM_WORLD;
    fclaw_mpi_init (&argc, &argv, mpicomm, lp);

    /* ------------------------------------
       Option handling
       ------------------------------------ */
    options = sc_options_new(argv[0]);

    /* Register [Options] and [clawpack46].  Very basic default values are set */
    fclaw_options_register(options,gparms);
    clawpack46_register_options(options,clawpack_parms);

    /* Read values first from fclaw2d_options.ini and then from the command line */
    retval = fclaw_options_read_from_file(options, lp);
    retval = retval || fclaw_options_parse_command_line (options,argc, argv, lp);

    /* post-process array options */
    fclaw_options_postprocess(gparms);
    clawpack46_postprocess_parms(clawpack_parms);

    /* Check final state of parameters.  Return from help message, if necessary. */
    retval = retval || fclaw_options_check (options, gparms, lp);
    retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms,lp);

    if (!retval)
    {
        /* Options are all okay; now run the example */
        fclaw_options_print_summary(options,lp);

        if (gparms->trapfpe == 1)
        {
            printf("Enabling floating point traps\n");
            feenableexcept(FE_INVALID);
        }

        /* ---------------------------------------------------------------
           Domain geometry
           -------------------------------------------------------------- */

        /* Map unit square to disk using mapc2m_disk.f */
        gparms->manifold = 0;
        conn = p4est_connectivity_new_unitsquare();
        cont = fclaw2d_map_new_nomap();

        domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

        if (gparms->verbosity > 0)
        {
            fclaw2d_domain_list_levels(domain, lp);
            fclaw2d_domain_list_neighbors(domain, lp);
        }

        /* ---------------------------------------------------------------
           Set domain data.
           --------------------------------------------------------------- */
        init_domain_data(domain);

        set_domain_parms(domain,gparms);
        set_clawpack_parms(domain,clawpack_parms);

        /* ---------------------------------------------------------------
           Define the solver and link in other problem/user specific
           routines
           --------------------------------------------------------------- */

        link_problem_setup(domain,swirl_problem_setup);

        swirl_link_solvers(domain);

        link_regrid_functions(domain,swirl_patch_tag4refinement,swirl_patch_tag4coarsening);

        /* ---------------------------------------------------------------
           Run
           --------------------------------------------------------------- */
        amrinit(&domain);
        amrrun(&domain);
        amrreset(&domain);

        /* This has to be in this scope */
        fclaw2d_map_destroy(cont);
    }
    else
    {
        /* Options are not okay;  Report error message and exit */
    }

    sc_options_destroy(options);
    fclaw_options_destroy_arrays(gparms);
    fclaw2d_clawpack_parms_delete(clawpack_parms);

    fclaw_mpi_finalize ();

    return 0;
}
