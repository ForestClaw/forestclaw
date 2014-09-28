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

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw2d_map.h>

#include "metric_user.H"

int
main (int argc, char **argv)
{
  int			lp;
  int                   example;
  MPI_Comm		mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;

#ifdef TRAPFPE
  printf("Enabling floating point traps\n");
  feenableexcept(FE_INVALID);
#endif

  lp = SC_LP_PRODUCTION;
  mpicomm = MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

#ifdef MPI_DEBUG
  /* This has to come after MPI has been initialized */
  fclaw2d_mpi_debug();
#endif

  /* ---------------------------------------------------------------
     Read in parameters and options
     --------------------------------------------------------------- */
  options = sc_options_new (argv[0]);

  sc_options_add_int (options, 0, "example", &example, 0,
                      "1 for Cartesian, " \
                      "2 for five patch square, " \
                      "3 for squared disk, " \
                      "4 for pillow disk, " \
                      "5 for pillow/five patch disk, " \
                      "6 for pillow sphere, " \
                      "7 for cubed sphere, "\
                      "8 for torus");

  /* Read parameters from .ini file */
  gparms = amr_options_new (options); // Sets default values

  /* Parse command line */
  amr_options_parse (options, argc, argv, lp);  // Reads options from a file

  /* Postprocess arrays */
  amr_postprocess_parms(gparms);
  amr_checkparms(gparms);

  /* ---------------------------------------------------------------
     Domain geometry
     --------------------------------------------------------------- */
  double alpha = 0.4;
  double scale[3];
  double shift[3];
  double rotate[2];
  set_default_transform(scale,shift,rotate);

  switch (example) {
  case 1:
      /* Map [0,1]x[0,1] to [-1,1],[-1,1] */
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_cart (scale, shift,rotate);
      break;
  case 2:
      /* Map [0,1],[0,1] to five patch square */
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_fivepatch (scale,shift,rotate, alpha);
      break;
  case 3:
      /* Map [0,1],[0,1] to squared disk */
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_squareddisk (scale,shift,rotate,alpha);
      break;
  case 4:
      /* Map [0,1],[0,1] to pillow disk */
      conn = p4est_connectivity_new_unitsquare ();
      cont = fclaw2d_map_new_pillowdisk (scale,shift,rotate);
      break;
  case 5:
      /* Map [0,1]x[0,1] to five patch --> pillow disk */
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_pillowdisk5 (scale,shift,rotate,alpha);
      break;
  case 6:
      /* Map [0,1]x[0,1] to five patch --> pillow disk */
      conn = p4est_connectivity_new_pillow ();
      cont = fclaw2d_map_new_pillowsphere (scale,shift,rotate);
      break;
  case 7:
      /* Map [0,1]x[0,1] to five patch --> pillow disk */
      conn = p4est_connectivity_new_cubed ();
      cont = fclaw2d_map_new_cubedsphere (scale,shift,rotate);
      break;
  case 8:
      /* Map [0,1]x[0,1] to five patch --> pillow disk */
      conn = p4est_connectivity_new_periodic ();
      cont = fclaw2d_map_new_torus (scale,shift,rotate,alpha);
      break;
  default:
      sc_abort_collective ("Parameter example must be 1 or 2");
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

  /* ----------------------------------------------------------
     to retrieve the context.  Note that this is only be used for
     passing the context to a C/C++ routine.  Do not expect to be
     able to access fields of the cont structure.
     ---------------------------------------------------------- */
  SET_CONTEXT(&cont);

  /* ---------------------------------------------------------- */

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

  /* Link functions that are called for the whole domain */
  link_problem_setup(domain,metric_setprob);
  link_run_diagnostics(domain,metric_diagnostics);

  /* Link other routines that need to be included. */
  metric_link_patch(domain);

  /* --------------------------------------------------
     Initialize and run the simulation
     -------------------------------------------------- */

  amrinit(&domain);

  if (gparms->run_diagnostics)
  {
      run_diagnostics(domain);
  }

  int iframe = 0;
  amrout(domain,iframe);


  /* amrrun(&domain); */

  amrreset(&domain);

  /* --------------------------------------------------
     Clean up.
     -------------------------------------------------- */
  amr_options_destroy(gparms);
  sc_options_destroy(options);

  fclaw_mpi_finalize ();

  return 0;
}
