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

#include "amr_single_step.h"
#include "fclaw2d_clawpack.H"

#include <amr_forestclaw.H>
#include <amr_utils.H>

#include "hemisphere_user.H"

int
main (int argc, char **argv)
{
  int			lp;
  MPI_Comm		mpicomm;
  sc_options_t          *options;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t* clawpack_parms;

  lp = SC_LP_PRODUCTION;
  mpicomm = MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* propose option handling as present in p4est/libsc */
  /* the option values live in amr_options, see amr_options.h */
  options = sc_options_new (argv[0]);

  /* Read in values from default .ini files */
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

  // ---------------------------------------------------------------
  // Domain geometry
  // ---------------------------------------------------------------
  domain = fclaw2d_domain_new_unitsquare(mpicomm,gparms->minlevel);

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

  sc_options_destroy (options);
  amr_options_destroy(gparms);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
