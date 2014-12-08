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
#include <fclaw_options.h>

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
    fclaw2d_map_data_t        smap_data, *map_data=&smap_data;

    int example, retval;
    const char* frames_string, *density_string;
    int *frames;
    double *density;
    double test_fpe;

    lp = SC_LP_PRODUCTION;
    mpicomm = sc_MPI_COMM_WORLD;
    fclaw_mpi_init (&argc, &argv, mpicomm, lp);

    /* Dictionary to store options */
    options = sc_options_new (argv[0]);

    /* [main] Register example option */
    sc_options_add_int (options, 0, "main:example", &example, 1,
                        "[main] Example = 1 or 2 [1]");

    fclaw_options_add_int_array(options,0,"main:frames",&frames_string,
                                "1 2", &frames,2,
                                "[main] Frames these frames [1 2]");

    fclaw_options_add_double_array(options,0,"main:density",&density_string,
                                   "1000.01 900.45",&density, 2,
                                   "[main] Density [1000.01 900.45]");

    /* [Options] General ForestClaw options */
    fclaw_options_register(options,gparms);

    /* [mapping] General mapping options (mi,mj,scale,shift,phi,theta) */
    fclaw2d_register_map_data(options,map_data); /* sets default values */

    /* [clawpack46] Options from Clawpack solver */
    clawpack46_register_options(options,clawpack_parms);

    /* Read from fclaw_options.ini */
    fclaw_options_read_from_file(options,lp);

    /* Override any values with command line values or from --inifile */
    retval = fclaw_options_parse_command_line(options,argc,argv,lp);

    /* convert array inputs */
    fclaw_options_postprocess(gparms);
    clawpack46_postprocess_parms(clawpack_parms);
    fclaw2d_options_postprocess_map_data(map_data);

    /* Post-process arrays created in [main] */
    fclaw_options_convert_int_array(frames_string,&frames,2);
    fclaw_options_convert_double_array(density_string,&density,2);

    /* Final check on options */
    retval = retval || fclaw_options_check(options,gparms,lp);
    retval = retval || test_parms_checkparms(example,lp);
    retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms,lp);

    if (!retval)
    {
        fclaw_options_print_summary(options,lp);

        if (gparms->mpi_debug == 1)
        {
            /* This doesn't work yet, since I don't have FCLAW_ENABLE_MPI defined */
            fclaw2d_mpi_debug();
        }

        if (gparms->trapfpe == 1)
        {
            feenableexcept(FE_INVALID);
        }

        printf("\n");
        switch (example)
        {
        case 1:
            printf("Running example 1 ... Done!\n");
            break;
        case 2:
            printf("Running example 2 ... \n");
            for(int i = 0; i < 2; i++)
            {
                printf("Density in frame %d : %f\n",frames[i],density[i]);
            }
            printf("Now let's trap a floating point error (this is not a bug!)\n");
            test_fpe = sqrt(-1.0);
            printf("Done!\n");
            break;
        default:
            SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
        }
        printf("\n");
    }
    fclaw_options_destroy_array((void*) frames);
    fclaw_options_destroy_array((void*) density);
    fclaw2d_map_destroy_arrays(map_data);
    fclaw2d_clawpack_parms_delete(clawpack_parms);
    fclaw_options_destroy_arrays(gparms);
    sc_options_destroy (options);

    fclaw_mpi_finalize ();

    return 0;
}
