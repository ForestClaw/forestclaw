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
  amr_options_t         samr_options, *gparms = &samr_options;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  lp = SC_LP_PRODUCTION;

  sc_init (mpicomm, 0, 0, NULL, lp);
  p4est_init (NULL, lp);

  /* propose option handling as present in p4est/libsc */
  /* the option values live in amr_options, see amr_options.h */
  options = sc_options_new (argv[0]);
  amr_options_register (options, gparms);
  amr_options_parse (options, gparms, argc, argv, lp);

  /* -----------------------------------------------------------------*/
  // This is set here just until we can read arrays in the SC library */
  Parser P;
  P.define(argc, argv);

  // This reads vector values into amr_options
  parse_ini_file(gparms);

  for (int j = 0; j < 2*SpaceDim; j++)
  {
      printf("mthbc[%d]= %d\n",j, gparms->mthbc[j]);
  }
  for (int j = 0; j < SpaceDim; j++)
  {
      printf("order[%d]= %d\n",j, gparms->order[j]);
  }
  for (int j = 0; j < gparms->mwaves; j++)
  {
      printf("mthlim[%d]= %d\n",j, gparms->mthlim[j]);
  }
  /* -----------------------------------------------------------------*/

  /* -----------------------------------------------------------------*/
  /* Sample user defined options */
  /* -----------------------------------------------------------------*/


  // Put this here so that we can read in the minimum level.
  // global_parms *gparms = new global_parms();
  // gparms->get_inputParams();

  int minlevel = gparms->minlevel;
  domain = fclaw2d_domain_new_unitsquare (mpicomm,minlevel);

  fclaw2d_domain_list_levels(domain, lp);
  fclaw2d_domain_list_neighbors(domain, lp);
  printf("\n\n");

  amrinit(&domain, gparms);
  amrrun(&domain);
  amrreset(&domain);

  // delete gparms;
  amr_options_delete (gparms);  // This doesn't do anything now...
  sc_options_destroy (options);
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
