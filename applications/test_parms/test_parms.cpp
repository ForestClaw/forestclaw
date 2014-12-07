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
#include <fclaw2d_options.h>

static int
    test_parms_checkparms (int example, int lp)
{
    if (example < 1 || example > 2) {
        fclaw2d_global_log (lp, "Option --example must be 1 or 2\n");
        return -1;
    }
    return 0;
}


int main (int argc, char **argv)
{
    int		              lp;
    sc_MPI_Comm               mpicomm;
    sc_options_t             *options;
    fclaw2d_clawpack_parms_t  sclawparms,  *clawpack_parms = &sclawparms;
    amr_options_t             sparms, *gparms = &sparms;

    int example, retval;
    int trapfpe, mpi_debug;
    double test_fpe;
    fclaw2d_map_data_t smap_data, *map_data=&smap_data;

    lp = SC_LP_PRODUCTION;
    mpicomm = sc_MPI_COMM_WORLD;
    fclaw_mpi_init (&argc, &argv, mpicomm, lp);

    /* Dictionary to store options */
    options = sc_options_new (argv[0]);

    /* [Example] Register example option */
    sc_options_add_int (options, 0, "main:example", &example, 1,
                        "[main] Example = 1 or 2 [1]");

    sc_options_add_bool (options, 0, "main:trapfpe", &trapfpe, 1,
                        "[main] Trap floating point exceptions [1]");

    sc_options_add_bool (options, 0, "main:mpi_debug", &mpi_debug, 0,
                        "[main] Start MPI debug process [0]");

    /* [mapping] General mapping parameters (mi,mj,scale,shift,phi,theta) */
    fclaw2d_register_map_data(options,map_data); /* sets default values */

    /* [forestclaw] General ForestClaw parameters */
    fclaw2d_register_options(options,gparms);

    /* [clawpack46] Parameters from Clawpack solver */
    clawpack46_register_options(options,clawpack_parms);

    fclaw_options_read_from_file(options,lp);  /* Read from fclaw_defaults.ini */

    /* Override any values with command line values or from --inifile */
    retval = fclaw_options_parse_command_line(options,argc,argv,lp);

    fclaw2d_postprocess_parms(gparms);  /* convert array inputs */
    clawpack46_postprocess_parms(clawpack_parms);  /* convert array inputs */

    retval = retval || fclaw2d_checkparms(options,gparms,lp);
    retval = retval || test_parms_checkparms(example,lp);
    retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms,lp);

    if (!retval)
    {
        fclaw_options_print_summary(options,lp);

        if (mpi_debug == 1)
        {
            /* This doesn't work yet, since I don't have FCLAW_ENABLE_MPI defined */
            fclaw2d_mpi_debug();
        }

        if (trapfpe == 1)
        {
            feenableexcept(FE_INVALID);
            test_fpe = sqrt(-1.0);
        }

        printf("\n");
        switch (example)
        {
        case 1:
            printf("Running example 1 ... Done!\n");
            break;
        case 2:
            printf("Running example 2 ... Done!\n");
            break;
        default:
            SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
        }
        printf("\n");
    }

    fclaw2d_clawpack_parms_delete(clawpack_parms);
    fclaw2d_options_destroy_arrays(gparms);
    sc_options_destroy (options);

    fclaw_mpi_finalize ();

    return 0;
}
