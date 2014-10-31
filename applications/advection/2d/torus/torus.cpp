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

#include <fclaw2d_clawpack.H>
#include <fclaw2d_map.h>
#include <p4est_connectivity.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw2d_map_query.h>

#include "torus_user.H"

static int
torus_checkparms (int example, int lp)
{
    if (example < 1 || example > 3) {
        fclaw2d_global_log (lp, "Option --example must be 1 or 2\n");
        return -1;
    }

    return 0;
}

int
main (int argc, char **argv)
{
  int                   lp;
  int                   retval;
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t* clawpack_parms;

  double pi, theta, phi;
  double rotate[2];
  double scale[3];
  double shift[3];
  double alpha;
  int example;
  int mi, mj, a,b;
  double R_earth;
  double lat[2];
  double longitude[2];


#ifdef TRAPFPE
  /* TODO: This should be done in fclaw_mpi_init? */
  printf("Enabling floating point traps\n");
  feenableexcept(FE_INVALID);
#endif

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

#ifdef MPI_DEBUG
  /* TODO: What is this for? Should be absorbed into fclaw_mpi_init. */
  /* This is code that suspends processing (with an infinite loop) so that
     multiple processes can be attached and debugged with GDB.  I am happy to
     move it.  Didn't know where... */
  fclaw2d_mpi_debug();
#endif
  /* ---------------------------------------------------------------
     Read parameters from .ini file, parse command line, and
     do parameter checking.
     -------------------------------------------------------------- */
  options = sc_options_new (argv[0]);

  sc_options_add_int (options, 0, "example", &example, 0,
                      "1=flat; 2 = torus");

  sc_options_add_double (options, 0, "theta", &theta, 0,
                         "Rotation angle theta (degrees) about z axis [0]");

  sc_options_add_double (options, 0, "phi", &phi, 0,
                         "Rotation angle phi (degrees) about x axis [0]");

  gparms = amr_options_new (options);
  clawpack_parms = fclaw2d_clawpack_parms_new(options);

  /* TODO: would it make sense to provide this code up to and excluding
   *       torus_checkparms into a reusable function? */
  /* parse command line; these values take precedence */
  amr_options_parse (options,argc, argv, lp);

  /* Any arrays are converted here */
  amr_postprocess_parms(gparms);
  fclaw2d_clawpack_postprocess_parms(clawpack_parms);

  /* Check final state of parameters.
   * The call to amr_checkparms2 also checks for a --help message.
   * We are exploiting C short-circuit boolean evaluation.
   */
  retval = amr_checkparms2 (options, gparms, lp);
  /* TODO: Convert this similarly to amr_checkparms2 */
  fclaw2d_clawpack_checkparms(clawpack_parms,gparms);
  /* TODO: example, phi, theta should live in a struct that can be passed */
  retval = retval || torus_checkparms (example, lp);
  if (!retval) {
     /* the do-the-work block. TODO: put everything below into a function */

  /* ---------------------------------------------------------------
     Domain geometry
     --------------------------------------------------------------- */
  pi = M_PI;
  set_default_transform(scale,shift,rotate);

  rotate[0] = pi*theta/180.0;
  rotate[1] = pi*phi/180.0;

  /* Ratio of inner radius to outer radius */
  alpha = 0.4;
  a = 1;  /* Periodicity */
  b = 1;

  /* default torus has outer radius 1 and inner radius alpha */
  switch (example)
  {
  case 1:
      a = 1;
      b = 0;
      mi = 4;
      mj = 4;
      scale[0] = (double) 2.0/mi;
      scale[1] = (double) 2.0/mj;
      shift[0] = -1;   /* [0,2]x[0,2] --> [-1,1]x[-1,1] */
      shift[1] = -1;
      conn = p4est_connectivity_new_brick(mi,mj,a,b);
      cont = fclaw2d_map_new_brick2(domain,mi,mj);
      break;
  case 2:
      /* Generally, we should have mj \approx \alpha*mi */
      mi = 2;
      mj = 1;
      conn = p4est_connectivity_new_brick(mi,mj,a,b);
      cont = fclaw2d_map_new_torus(scale,shift,rotate,alpha,mi,mj);
      break;
  case 3:
      /* Lat-long example */
      mi = 2;
      mj = 1;
      a = 1;
      b = 0;
      longitude[0] = 0; /* x-coordinate */
      longitude[1] = 360;  /* if a == 1, long[1] will be computed as [0] + 180 */
      lat[0] = -80;  /* y-coordinate */
      lat[1] = 80;
      R_earth = 6378.1;
      scale[0] = 1;
      scale[1] = 1;
      scale[2] = 1;
      conn = p4est_connectivity_new_brick(mi,mj,a,b);
      cont = fclaw2d_map_new_latlong(scale,shift,rotate,lat,longitude,mi,mj,a,b);
      break;
  default:
      SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);


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
  set_clawpack_parms(domain,clawpack_parms);

  link_problem_setup(domain,torus_problem_setup);

  torus_link_solvers(domain);

  link_regrid_functions(domain,
                        torus_patch_tag4refinement,
                        torus_patch_tag4coarsening);


  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  } /* this bracket ends the do_the_work block */

  sc_options_destroy (options);         /* this could be moved up */
  amr_options_destroy (gparms);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
