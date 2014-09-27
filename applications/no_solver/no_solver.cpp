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

// #include <amr_single_step.h>
// #include <fclaw2d_clawpack.H>
// #include <amr_forestclaw.H>
// #include <amr_utils.H>

#include <fclaw2d_map.h>
#include <p4est_connectivity.h>
#include <fclaw2d_map_query.h>


#include "no_solver_user.H"

int
main (int argc, char **argv)
{
  int		        lp;
  int                   example;
  sc_MPI_Comm           mpicomm;
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
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

#ifdef MPI_DEBUG
  /* This has to come after MPI has been initialized */
  fclaw2d_mpi_debug();
#endif


  /* ---------------------------------------------------------------
     Read parameters from .ini file, parse command line, and
     do parameter checking.
     -------------------------------------------------------------- */
  options = sc_options_new(argv[0]);

  sc_options_add_int (options, 0, "example", &example, 0,
                      "0 no mapping (use [ax,bx]x[ay,by], " \
                      "1 for Cartesian," \
                      "2 for five patch square");

  /* Register default parameters and any solver parameters */
  gparms = amr_options_new(options);

  amr_options_parse(options,argc,argv,lp);

  /* Any arrays are converted here */
  amr_postprocess_parms(gparms);

  /* Check final state of parameters */
  amr_checkparms(gparms);

  /* ---------------------------------------------------------------
     Domain geometry
     -------------------------------------------------------------- */

  double alpha = 0.5;

  double scale[3];
  double shift[3];
  double rotate[2];
  set_default_transform(scale,shift,rotate);

  switch (example) {
  case 0:
      /* Don't use a mapping.  [ax,ay]x[ay,by] will be used instead */
      gparms->manifold = 0;
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_nomap ();
      break;
  case 1:
      /* Map unit square to disk using mapc2m_disk.f */
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_cart (scale, shift, rotate);
      break;
  case 2:
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_fivepatch (scale,shift,rotate,alpha);
      break;
  case 3:
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_pillowdisk (scale,shift,rotate);
      break;
  case 4:
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_pillowfivepatch (scale,shift,rotate,alpha);
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
  set_domain_parms(domain,gparms);

  /* ---------------------------------------------------------------
     Don't do any solve or run.  This is just to test that we can
     compile without any solver routines.
     --------------------------------------------------------------- */

  no_solver_linker(domain);

  /* ---------------------------------------------------------------
     Initialize and run (but with out updating anything)
     --------------------------------------------------------------- */
  amrinit(&domain);
  amrrun(&domain);  /* Nothing should happen */
  amrreset(&domain);

  /* ---------------------------------------------------------------
     Clean up
     --------------------------------------------------------------- */

  sc_options_destroy(options);         /* this could be moved up */
  amr_options_destroy(gparms);

  fclaw_mpi_finalize ();

  return 0;
}
