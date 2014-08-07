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

#include "cubed_sphere_user.H"

int
main (int argc, char **argv)
{
  int			lp;
  int                   example;
  double                cubed_R1, cubed_R2;
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL;
  fclaw2d_domain_t	*domain = NULL;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t* clawpack_parms;

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* ---------------------------------------------------------------
     Read in parameters and options
     --------------------------------------------------------------- */
  options = sc_options_new (argv[0]);
  sc_options_add_int (options, 0, "example", &example, 0,
                      "2 for disk, 3 for sphere");
  sc_options_add_double (options, 0, "cubed_R1", &cubed_R1, -1.,
                         "outer radius of cubed disk and cubed sphere");
  sc_options_add_double (options, 0, "cubed_R2", &cubed_R2, -1.,
                         "circumcircle radius of inner cube of cubed disk");

  /* Read parameters from .ini file */
  gparms = amr_options_new (options); // Sets default values
  clawpack_parms = fclaw2d_clawpack_parms_new(options);

  amr_options_parse (options, argc, argv, lp);  // Reads options from a file

  amr_postprocess_parms (gparms);
  fclaw2d_clawpack_postprocess_parms(clawpack_parms);

  /* Read in clawpack parms, process and check them */
  amr_checkparms (gparms);
  fclaw2d_clawpack_checkparms(clawpack_parms,gparms);

  /* ---------------------------------------------------------------
     Domain geometry
     --------------------------------------------------------------- */

  // Queries for sphere grids.
  int query_results[FCLAW2D_MAP_QUERY_LAST];
  query_results[0] = 1;
  query_results[1] = 0;
  query_results[2] = 0;
  query_results[3] = 1;
  query_results[4] = 0;
  query_results[5] = 0;
  query_results[6] = 0;


  /*
    // Save disk for another example.
  case 2:
      if (!(0. < cubed_R2 && cubed_R2 < cubed_R1)) {
          sc_abort_collective
          ("Parameters 0 < cubed_R2 < cubed_R1 required for cubed disk");
      }
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_disk (cubed_R1, cubed_R2);
      break;
  */

  switch (example) {
  case 1:
      set_maptype_pillowsphere_();
      conn = p4est_connectivity_new_pillow();
      cont = fclaw2d_map_new_fortran(mapc2m_,query_results);
      break;
  case 2:
      if (!(0. < cubed_R1)) {
        sc_abort_collective ("Parameter 0 < cubed_R1 required for cubed sphere");
      }
      conn = p4est_connectivity_new_cubed ();
      cont = fclaw2d_map_new_csphere (cubed_R1);
      break;
  case 3:
      set_maptype_cubedsphere_();
      conn = p4est_connectivity_new_cubed();
      cont = fclaw2d_map_new_fortran(mapc2m_,query_results);
      break;
    default:
      sc_abort_collective ("Parameter example must be 1, 3 or 4");
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
  link_problem_setup(domain,cubed_sphere_problem_setup);

  cubed_sphere_link_solvers(domain);

  /* --------------------------------------------------
     Initialize and run the simulation
     -------------------------------------------------- */
  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  /* --------------------------------------------------
     Clean up.
     -------------------------------------------------- */
  fclaw2d_map_destroy (cont);
  amr_options_destroy(gparms);
  sc_options_destroy(options);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
