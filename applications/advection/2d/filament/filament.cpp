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

#include "filament_user.H"

int
main (int argc, char **argv)
{
  int		        lp;
  int                   example;
  double                R, xcenter,ycenter;  // square has width w=2*R
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t  *clawpack_parms;

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* ---------------------------------------------------------------
     Read parameters from .ini file, parse command line, and
     do parameter checking.
     -------------------------------------------------------------- */
  options = sc_options_new(argv[0]);
  sc_options_add_int (options, 0, "example", &example, 0,
                      "1 for pillow 'disk' (shifted and scaled), "\
                      "2 for 5-square disk (scaled), " \
                      "3 for unit square (scaled)");

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
      /* Setting up a mapping using an existing mapc2m.f function */
      set_maptype_disk_();
      int query_results[FCLAW2D_MAP_QUERY_LAST];
      query_results[0] = 1;  // IS_USED
      query_results[1] = 0;  // IS_SCALESHIFT
      query_results[2] = 0;  // IS_AFFINE
      query_results[3] = 1;  // IS_NONLINEAR
      query_results[4] = 0;  // IS_GRAPH
      query_results[5] = 1;  // IS_PLANAR
      query_results[6] = 0;  // IS_ALIGNED
      query_results[7] = 1;  // IS_FLAT
      query_results[8] = 0;  // IS_PILLOWDISK

      /* Map unit square to disk using mapc2m_disk.f */
      conn = p4est_connectivity_new_unitsquare();
      cont = fclaw2d_map_new_fortran(mapc2m_disk,query_results);
      break;
  case 2:
      if (!(0. < R)) {
        sc_abort_collective ("Parameter 0 < R required for 5-square disk");
      }
      conn = p4est_connectivity_new_disk ();
      cont = fclaw2d_map_new_disk (R);
      break;
  case 3:
      conn = p4est_connectivity_new_pillowdisk();
      cont = fclaw2d_map_new_pillowdisk();
      break;
    default:
      sc_abort_collective ("Parameter example must be 1, 3 or 4");
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);




  domain = fclaw2d_domain_new_unitsquare (mpicomm, gparms->minlevel);

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

  swirl_link_solvers(domain);

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
