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

#include "circles_user.H"

int
main (int argc, char **argv)
{
  int			lp;
  MPI_Comm		mpicomm;
  fclaw2d_domain_t	*domain;
  fclaw_options_t         samr_options, *gparms = &samr_options;

  lp = SC_LP_PRODUCTION;
  mpicomm = MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* ---------------------------------------------------------------
     Read in parameters and options
     --------------------------------------------------------------- */
  options = sc_options_new (argv[0]);

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
  /* For sphere */
  // domain = fclaw2d_domain_new_twosphere (mpicomm,gparms->minlevel);
  domain = fclaw2d_domain_new_unitsquare (mpicomm, gparms->minlevel);

  fclaw2d_domain_list_levels(domain, lp);
  fclaw2d_domain_list_neighbors(domain, lp);

  /* ---------------------------------------------------------------
     Set domain data.
     --------------------------------------------------------------- */
  init_domain_data(domain);

  /* Store parameters */
  set_domain_parms(domain,gparms);

  link_problem_setup(domain,sphere_setprob);
  link_regrid_functions(domain,sphere_patch_tag4refinement,
                        sphere_patch_tag4coarsening);
  link_output_functions(domain,metric_write_header, metric_write_output);

  // Link other routines that need to be included.
  sphere_link_patch(domain);

  /* --------------------------------------------------
     Initialize and run the simulation
     -------------------------------------------------- */
  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  /* --------------------------------------------------
     Clean up.
     -------------------------------------------------- */
  amr_options_destroy(gparms);
  sc_options_destroy(options);

  fclaw_mpi_finalize ();

  return 0;
}
