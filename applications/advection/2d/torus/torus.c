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

#include <fclaw2d_capi.h>
#include <fclaw2d_convenience.h>
#include <amr_single_step.H>
#include <amr_waveprop.H>

int
main (int argc, char **argv)
{
  int                   lp;
  MPI_Comm              mpicomm;
  sc_options_t          *options;
  fclaw2d_domain_t	*domain;
  amr_options_t         *gparms;

  lp = SC_LP_PRODUCTION;
  mpicomm = MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  /* the option values live in amr_options_t * gparms; see amr_options.h */
  options = sc_options_new (argv[0]);
  /* set default options and read default option file */
  gparms = amr_options_new (options);
  /* parse command line; these values take precedence */
  amr_options_parse (options, gparms, argc, argv, lp);

  /* ---------------------------------------------------------------
     Domain geometry
     --------------------------------------------------------------- */
  domain = fclaw2d_domain_new_unitsquare (mpicomm, gparms->minlevel);

  fclaw2d_domain_list_levels(domain, lp);
  fclaw2d_domain_list_neighbors(domain, lp);

  /* ---------------------------------------------------------------
     Set domain data.
     --------------------------------------------------------------- */
  fclaw2d_allocate_domain_data (domain, gparms,
                                &fclaw2d_single_step,
                                &fclaw2d_waveprop_update);

  amrinit(&domain);
  amrrun(&domain);
  amrreset(&domain);

  sc_options_destroy (options);         /* this could be moved up */
  amr_options_destroy (gparms);

  fclaw_mpi_finalize ();

  return 0;
}
