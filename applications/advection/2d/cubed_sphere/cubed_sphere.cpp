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

#include "cubed_sphere_user.H"

int
main (int argc, char **argv)
{
  int			lp;
  int                   example;
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL;
  fclaw2d_domain_t	*domain = NULL;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t* clawpack_parms;

  double scale, theta, phi;
  double rotate[2];

#ifdef TRAPFPE
  printf("Enabling floating point traps\n");
  feenableexcept(FE_INVALID);
#endif

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* ---------------------------------------------------------------
     Read in parameters and options
     --------------------------------------------------------------- */
  options = sc_options_new (argv[0]);
  sc_options_add_int (options, 0, "example", &example, 0,
                      "1 for pillow grid, "\
                      "2 for cubed sphere ");

  sc_options_add_double (options, 0, "scale", &scale, 1.0,
                         "Scale unit sphere (e.g. set radius [1])");
  sc_options_add_double (options, 0, "theta", &theta, 0,
                         "Rotation angle theta (degrees) about z axis [0]");

  sc_options_add_double (options, 0, "phi", &phi, 0,
                         "Rotation angle phi (degrees) about x axis [0]");

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

  double pi = M_PI;
  rotate[0] = pi*theta/180.0;
  rotate[1] = pi*phi/180.0;

  switch (example) {
  case 1:
      conn = p4est_connectivity_new_pillow();
      cont = fclaw2d_map_new_pillowsphere(rotate,scale);
      break;
  case 2:
      conn = p4est_connectivity_new_cubed();
      cont = fclaw2d_map_new_cubedsphere(rotate,scale);
      break;
    default:
      sc_abort_collective ("Parameter example must be 1 (pillow sphere) or 2 (cubed sphere)");
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

  /* ----------------------------------------------------------
     Set mapping context for Fortran.  The context will now be
     available via get_context(), as a integer*8.  This argument
     will show up as a (fclaw2d_map_context_t**) in C/C++.
     From Fortran, use

     c      # From fortran :
            integer*8 cont
            integer blockno

     c      # .......

            cont = get_context()
            blockno = get_block()

     c      # Call the mapping function
            call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)

     c      # .......

     to retrieve the context.  Note that this is only be used for
     passing the context to a C/C++ routine.  Do not expect to be
     able to access fields of the cont structure.
     ---------------------------------------------------------- */
    SET_CONTEXT(&cont);

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
