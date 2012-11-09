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

#include "fclaw2d_convenience.h"
#include "amr_options.h"
#include "amr_forestclaw.H"

#include "parser.H"

// This needs to go away.  The p4est namespace should not be used directly.
#include <p4est.h>


int
main (int argc, char **argv)
{
  int			mpiret;
  int			lp;
  MPI_Comm		mpicomm;
  sc_options_t          *options;
  fclaw2d_domain_t	*domain;
  amr_options_t         *gparms;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  lp = SC_LP_PRODUCTION;

  sc_init (mpicomm, 0, 0, NULL, lp);
  p4est_init (NULL, lp);

  /* propose option handling as present in p4est/libsc */
  /* the option values live in amr_options, see amr_options.h */
  options = sc_options_new (argv[0]);
  gparms = amr_options_new (options);
  amr_options_parse (options, gparms, argc, argv, lp);

  /* -----------------------------------------------------------------*/
#if 0
  // This is set here just until we can read arrays in the SC library */
  Parser P;
  P.define(argc, argv);

  // This reads vector values into amr_options
  parse_ini_options(gparms);
#endif

  for (int j = 0; j < SpaceDim; j++)
  {
      printf("order[%d]= %d\n",j, gparms->order[j]);
  }

  exit (1);
  for (int j = 0; j < 2*SpaceDim; j++)
  {
      printf("mthbc[%d]= %d\n",j, gparms->mthbc[j]);
  }
  for (int j = 0; j < gparms->mwaves; j++)
  {
      printf("mthlim[%d]= %d\n",j, gparms->mthlim[j]);
  }

  parse_ini_output(gparms,gparms->outstyle);

  /*
  sc_options_print_usage(sc_package_id,lp,options,"");
  exit(1);
  */

  /* -----------------------------------------------------------------*/

  /* -----------------------------------------------------------------*/
  /* Sample user defined options */
  /* -----------------------------------------------------------------*/

  domain = fclaw2d_domain_new_unitsquare (mpicomm, gparms->minlevel);

  fclaw2d_domain_list_levels(domain, lp);
  fclaw2d_domain_list_neighbors(domain, lp);
  printf("\n\n");

  amrinit(&domain, gparms);
  amrrun(&domain);
  amrreset(&domain);

  amr_options_destroy (gparms);
  sc_options_destroy (options);
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
