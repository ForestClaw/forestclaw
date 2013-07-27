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

// This needs to go away.  The p4est namespace should not be used directly.
#include <p4est.h>

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "amr_manyclaw.H"
#include "amr_options.h"

#include "quadrants_user.H"

int
main (int argc, char **argv)
{
  int		        lp;
  MPI_Comm              mpicomm;
  sc_options_t          *options;
  fclaw2d_domain_t	*domain;
  amr_options_t         *gparms;
  amr_manyclaw_parms_t  *manyclaw_parms;

  lp = SC_LP_PRODUCTION;
  mpicomm = MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* ---------------------------------------------------------------
     Read parameters from .ini file, parse command line, and
     do parameter checking.
     -------------------------------------------------------------- */
  options = sc_options_new(argv[0]);

  /* Register default parameters and any solver parameters */
  gparms = amr_options_new(options);
  manyclaw_parms = amr_manyclaw_parms_new(options);

  /* Parse any command line arguments.  Argument gparms is no longer needed
     as an input since the array conversion is now done in 'postprocess' */

  amr_options_parse(options,argc,argv,lp);

  /* Any arrays are converted here */
  amr_postprocess_parms(gparms);
  amr_manyclaw_postprocess_parms(manyclaw_parms);

  /* Check final state of parameters */
  amr_checkparms(gparms);
  amr_manyclaw_checkparms(manyclaw_parms,gparms);

  /* ---------------------------------------------------------------
     Domain geometry
     -------------------------------------------------------------- */
  domain = fclaw2d_domain_new_unitsquare (mpicomm, gparms->minlevel);

  fclaw2d_domain_list_levels(domain, lp);
  fclaw2d_domain_list_neighbors(domain, lp);

  /* ---------------------------------------------------------------
     Set domain data.
     --------------------------------------------------------------- */
  init_domain_data(domain);

  set_domain_parms(domain,gparms);
  set_manyclaw_parms(domain,manyclaw_parms);

  /* ---------------------------------------------------------------
     Define the solver and link in other problem/user specific
     routines
     --------------------------------------------------------------- */

  /* Using user defined functions just to demonstrate how one might setup
     something that depends on more than one solver (although only one is used
     here) */
  link_problem_setup(domain,quadrants_problem_setup);

  quadrants_link_solvers(domain);

  link_regrid_functions(domain,quadrants_patch_tag4refinement,quadrants_patch_tag4coarsening);

  /* ---------------------------------------------------------------
     Run
     --------------------------------------------------------------- */

  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  sc_options_destroy(options);         /* this could be moved up */
  amr_options_destroy(gparms);
  amr_manyclaw_parms_delete(manyclaw_parms);

  fclaw_mpi_finalize ();

  return 0;
}
