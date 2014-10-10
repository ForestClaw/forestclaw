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


int
main (int argc, char **argv)
{
  int		        lp;
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);


  /* ----------------------------------------------------------
     Read in command line options
     ---------------------------------------------------------- */
  int example;
  options = sc_options_new (argv[0]);
  sc_options_add_int (options, 0, "example", &example, 0,
                      "1 for AMR options, " \
                      "2 for Clawpack parms");

  /* ---------------------------------------------------------------
     Read parameters from .ini file, parse command line, and
     do parameter checking.
     -------------------------------------------------------------- */
  amr_options_parse(options,argc,argv,lp);


  switch (example) {
  case 1:
      amr_options_t  *gparms;

      /* Register default parameters and any solver parameters */
      gparms = amr_options_new(options);
      amr_options_parse(options,argc,argv,lp);
      amr_postprocess_parms(gparms);

      break;
  case 2:
      fclaw2d_clawpack_parms_t  *clawpack_parms;

      clawpack_parms = fclaw2d_clawpack_parms_new(options);

      /* Command line arguments */
      amr_options_parse(options,argc,argv,lp);

      fclaw2d_clawpack_postprocess_parms(clawpack_parms);

      fclaw2d_clawpack_parms_delete(clawpack_parms);

      break;
  default:
      sc_abort_collective ("Use command line --example=1 or --example=2");
  }


  /* -----------------------------------------------
     Skip everything that actually does something...
     ----------------------------------------------- */

  sc_options_destroy (options);

  fclaw_mpi_finalize ();

  return 0;
}
