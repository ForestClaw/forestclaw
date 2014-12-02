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
#include <fclaw2d_map_query.h>

#include "hemisphere_user.H"

int
main (int argc, char **argv)
{
  int			lp;
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t* clawpack_parms;

  double theta;

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

#ifdef MPI_DEBUG
  /* this has to go after MPI has been initialized */
  fclaw2d_mpi_debug();
#endif


  /* ----------------------------------------------------------
     Read in command line options
     ---------------------------------------------------------- */
  int example;
  options = sc_options_new (argv[0]);
  sc_options_add_int (options, 0, "example", &example, 0,
                      "1 for Cartesian, " \
                      "2 for five patch square");

  sc_options_add_double (options, 0, "theta", &theta, 0,
                         "Rotation angle theta (degrees) about z axis [0]");

  /* ----------------------------------------------------------
     Read in values from .ini files.  These are overwritten by
     command line options read above.
     ---------------------------------------------------------- */
  gparms = amr_options_new (options);
  clawpack_parms = fclaw2d_clawpack_parms_new(options);

  /* Parse command line for any modifications */
  amr_options_parse (options, argc, argv, lp);  // Reads options from a file

  /* Postprocess array inputs */
  amr_postprocess_parms(gparms);
  fclaw2d_clawpack_postprocess_parms(clawpack_parms);

  /* Verify inputs */
  amr_checkparms(gparms);
  fclaw2d_clawpack_checkparms(clawpack_parms,gparms);

  /* ---------------------------------------------------------------
     Floating point traps
     -------------------------------------------------------------- */
  if (gparms->trapfpe == 1)
  {
      printf("Enabling floating point traps\n");
      feenableexcept(FE_INVALID);
  }

  /* ---------------------------------------------------------------
     Set up the domain geometry
     --------------------------------------------------------------- */
  double alpha = 0.5;
  double scale[3];
  double shift[3];
  double rotate[2];
  set_default_transform(scale,shift,rotate);
  rotate[0] = theta;
  rotate[1] = 0;


  switch (example) {
  case 1:
      /* Map unit square to disk using mapc2m_disk.f */
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_pillowsphere(scale,shift,rotate);
      break;
  case 2:
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_pillowsphere5(scale,shift,rotate,alpha);
      break;
  default:
      sc_abort_collective ("Parameter example must be 1 or 2");
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

  /* ---------------------------------------------------------------
     Print out domain info
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

  /* Store domain parameters */
  set_domain_parms(domain,gparms);
  set_clawpack_parms(domain,clawpack_parms);

/* ---------------------------------------------
   Define the solver
   --------------------------------------------- */

  /* Link waveprop solvers to domain */
  link_problem_setup(domain,hemisphere_problem_setup);

  hemisphere_link_solvers(domain);
  link_regrid_functions(domain,hemisphere_patch_tag4refinement,
                        hemisphere_patch_tag4coarsening);

  /* --------------------------------------------------
     Initialize and run the simulation
     -------------------------------------------------- */
  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  fclaw2d_map_destroy(cont);
  sc_options_destroy (options);
  amr_options_destroy(gparms);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
