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
#include <fclaw2d_clawpack.H>
#include <fclaw2d_map.h>
#include <p4est_connectivity.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw_options.h>

#include <fclaw2d_map_query.h>

#include "sphere_user.H"

static int
    sphere_checkparms (int example, int lp)
{
    if (example < 1 || example > 2) {
        fclaw2d_global_log (lp, "Option --example must be 1 or 2\n");
        return -1;
    }
    return 0;
}


int main (int argc, char **argv)
{
    int			    lp;
    sc_MPI_Comm             mpicomm;
    sc_options_t            *options;
    p4est_connectivity_t    *conn = NULL;
    fclaw2d_map_context_t   *cont = NULL;

    fclaw2d_domain_t	    *domain = NULL;

    amr_options_t            samr_options, *gparms = &samr_options;
    fclaw2d_clawpack_parms_t sclawpack_parms, *clawpack_parms = &sclawpack_parms;

    /* double theta, phi; */
    fclaw2d_map_data_t smap_data, *map_data = &smap_data;
    double pi = M_PI;
    double rotate[2];

    int example, retval;

    lp = SC_LP_PRODUCTION;
    mpicomm = sc_MPI_COMM_WORLD;
    fclaw_mpi_init (&argc, &argv, mpicomm, lp);

    /* ---------------------------------------------------------------
       Read in parameters and options
       --------------------------------------------------------------- */
    options = sc_options_new (argv[0]);

    /* [main] options */
    sc_options_add_int (options, 0, "main:example", &example, 0,
                        "[main] 1 for pillow grid, "\
                        "2 for cubed sphere ");

    /* [Options] General mapping functions */
    fclaw_options_register (options,gparms);

    /* [mapping] General mapping functions */
    fclaw2d_register_map_data(options,map_data); /* sets default values */

    /* [clawpack46] Clawpack solver options */
    clawpack46_register_options(options,clawpack_parms);

    retval = clawpack46_options_read_from_file(options,lp);
    retval = retval || fclaw_options_read_from_file(options,lp);

    retval = retval || fclaw_options_parse_command_line (options, argc, argv, lp);

    /* Post-process any array input */
    fclaw_options_postprocess (gparms);
    clawpack46_postprocess_parms(clawpack_parms);

    retval = retval || fclaw_options_check (options, gparms, lp);
    retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms,lp);
    retval = retval || sphere_checkparms (example, lp);

    if (!retval)
    {
        fclaw_options_print_summary(options,lp);

        if (gparms->trapfpe == 1)
        {
            printf("Enabling floating point traps\n");
            feenableexcept(FE_INVALID);
        }

        if (gparms->mpi_debug == 1)
        {
            fclaw2d_mpi_debug();
        }

        rotate[0] = pi*map_data->theta/180.0;
        rotate[1] = pi*map_data->phi/180.0;

        switch (example) {
        case 1:
            conn = p4est_connectivity_new_pillow();
            cont = fclaw2d_map_new_pillowsphere(map_data->scale,map_data->shift,rotate);
            break;
        case 2:
            conn = p4est_connectivity_new_cubed();
            cont = fclaw2d_map_new_cubedsphere(map_data->scale,map_data->shift,rotate);
            break;
        default:
            SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
        }

        domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

        /* ---------------------------------------------------------------
           Print some diagnostics.  TODO: Should this really be in main()?
           --------------------------------------------------------------- */

        if (gparms->verbosity > 0)
        {
            fclaw2d_domain_list_levels(domain, lp);
            fclaw2d_domain_list_neighbors(domain, lp);
        }

        /* ---------------------------------------------------------------
           Set domain data.
           --------------------------------------------------------------- */
        init_domain_data(domain);

        /* Store parameters */
        set_domain_parms(domain,gparms);
        set_clawpack_parms(domain,clawpack_parms);

        /* Link solvers to the domain */
        link_problem_setup(domain,sphere_problem_setup);

        sphere_link_solvers(domain);

        /* --------------------------------------------------
           Initialize and run the simulation
           -------------------------------------------------- */
        amrinit(&domain);
        amrrun(&domain);
        amrreset(&domain);

        /* --------------------------------------------------
           Clean up the mapping context.
           -------------------------------------------------- */
        fclaw2d_map_destroy (cont);
    }

    fclaw_options_destroy_arrays(gparms);  /* Delete any array options */
    sc_options_destroy(options);
    fclaw2d_clawpack_parms_delete(clawpack_parms);

    fclaw_mpi_finalize ();

    return 0;
}
