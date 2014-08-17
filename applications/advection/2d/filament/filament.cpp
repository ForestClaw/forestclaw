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

#include "filament_user.H"

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
  fclaw2d_clawpack_parms_t  *clawpack_parms;

  double rotate[2];
  double scale;
  double theta, phi;

#ifdef TRAPFPE
  printf("Enabling floating point traps\n");
  feenableexcept(FE_INVALID);
#endif


  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* ---------------------------------------------------------------
     Read parameters from .ini file, parse command line, and
     do parameter checking.
     -------------------------------------------------------------- */
  options = sc_options_new(argv[0]);
  sc_options_add_int (options, 0, "example", &example, 0,
                      "1 for pillow disk (scaled)," \
                      "2 for 5-square disk (scaled), " \
                      "3 for unit square (scaled)");

  sc_options_add_double (options, 0, "scale", &scale, 1.0,
                         "Scale unit sphere (e.g. set radius [1])");

  sc_options_add_double (options, 0, "theta", &theta, 0,
                         "Rotation angle theta (degrees) about z axis [0]");

  sc_options_add_double (options, 0, "phi", &phi, 0,
                         "Rotation angle phi (degrees) about x axis [0]");


  /* Register default parameters and any solver parameters */
  gparms = amr_options_new(options);
  clawpack_parms = fclaw2d_clawpack_parms_new(options);

  /* Parse any command line arguments.  Argument gparms is no longer needed
     as an input since the array conversion is now done in 'postprocess' */

  amr_options_parse(options,argc,argv,lp);

  /* Any arrays are converted here */
  amr_postprocess_parms(gparms);
  fclaw2d_clawpack_postprocess_parms(clawpack_parms);

  /* Check final state of parameters */
  amr_checkparms(gparms);
  fclaw2d_clawpack_checkparms(clawpack_parms,gparms);

  /* ---------------------------------------------------------------
     Domain geometry
     -------------------------------------------------------------- */
  double pi = 3.141592653589793;
  rotate[0] = pi*theta/180.0;
  rotate[1] = pi*phi/180.0;

  double R1 = 1.0;
  double R2 = 0.35;


  switch (example) {
  case 1:
      /* Map unit square to disk using mapc2m_disk.f */
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_pillowdisk(rotate,scale);
      break;
  case 2:
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_squareddisk (rotate,scale,R1,R2);
      break;
  case 3:
      /* in [0,1]x[0,1] */
      scale = 2.0;
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_cart(scale);
      break;
    default:
      sc_abort_collective ("Parameter example must be 1, 2 or 3");
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

  /* ----------------------------------------------------------
     to retrieve the context.  Note that this is only be used for
     passing the context to a C/C++ routine.  Do not expect to be
     able to access fields of the cont structure.
     ---------------------------------------------------------- */
    SET_CONTEXT(&cont);


  /* ----------------------------------------------------------
     ---------------------------------------------------------- */
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

  /* Using user defined functions just to demonstrate how one might setup
     something that depends on more than one solver (although only one is used
     here) */
  link_problem_setup(domain,filament_problem_setup);

  filament_link_solvers(domain);

  link_regrid_functions(domain,
                        filament_patch_tag4refinement,
                        filament_patch_tag4coarsening);

  /* ---------------------------------------------------------------
     Run
     --------------------------------------------------------------- */

  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  sc_options_destroy(options);         /* this could be moved up */
  amr_options_destroy(gparms);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
