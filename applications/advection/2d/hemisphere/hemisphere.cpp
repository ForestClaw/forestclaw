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

#include <fclaw2d_single_step.h>
#include <fclaw2d_waveprop.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>

int
main (int argc, char **argv)
{
  int			lp;
  MPI_Comm		mpicomm;
  sc_options_t          *options;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;

  lp = SC_LP_PRODUCTION;
  mpicomm = MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* propose option handling as present in p4est/libsc */
  /* the option values live in amr_options, see amr_options.h */
  options = sc_options_new (argv[0]);
  gparms = amr_options_new (options); // Sets default values
  amr_options_parse (options, gparms, argc, argv, lp);  // Reads options from a file

  // ---------------------------------------------------------------
  // Domain geometry
  // ---------------------------------------------------------------
  domain = fclaw2d_domain_new_unitsquare(mpicomm,gparms->minlevel);

  fclaw2d_domain_list_levels(domain, lp);
  fclaw2d_domain_list_neighbors(domain, lp);
  printf("\n\n");

  /* ---------------------------------------------------------------
     Set domain data.
     --------------------------------------------------------------- */
  allocate_user_data(domain);       // allocate all data for the domain.

  fclaw2d_domain_data_t *ddata = get_domain_data(domain);
  ddata->amropts = gparms;

/* ---------------------------------------------
   Define the solver
   --------------------------------------------- */
  // This simplified view of the world will have to change
  // when it comes time to think about solving, say,
  // advection-diffusion.
  ddata->f_level_advance = &fclaw2d_single_step;
  ddata->f_single_step_patch = &fclaw2d_waveprop_update;

  /* --------------------------------------------------
     Initialize and run the simulation
     -------------------------------------------------- */
  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  amr_options_destroy(gparms);
  sc_options_destroy (options);

  fclaw_mpi_finalize ();

  return 0;
}
