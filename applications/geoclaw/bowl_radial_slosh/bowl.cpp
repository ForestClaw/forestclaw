/*
  Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "radial/radial_user.h"
#include "slosh/slosh_user.h"

#include <fclaw_include_all.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fc2d_geoclaw.h>
#include <fc2d_geoclaw_options.h>

int
main (int argc, char **argv)
{
    /* Initialize application */
    fclaw_app_t *app = fclaw_app_new (&argc, &argv, NULL);

    /* Options */

    radial_user_options_t       *radial_user_opt;
    fclaw_options_t             *radial_gparms;
    fclaw_clawpatch_options_t *radial_clawpatchopt;
    fc2d_geoclaw_options_t      *radial_geoclawopt;

    slosh_user_options_t        *slosh_user_opt;
    fclaw_options_t             *slosh_fclaw_opt;
    fclaw_clawpatch_options_t *sloshclawpatch_opt;
    fc2d_geoclaw_options_t      *slosh_geo_opt;


    radial_gparms                   = fclaw_options_register(app, "radial",           "fclaw_options.ini");
    radial_clawpatchopt = fclaw_clawpatch_2d_options_register(app, "radial-clawpatch", "fclaw_options.ini");
    radial_geoclawopt        = fc2d_geoclaw_options_register(app, "radial-geoclaw",   "fclaw_options.ini");
    radial_user_opt =                radial_options_register(app, "radial-user",      "fclaw_options.ini");  

    slosh_fclaw_opt =                   fclaw_options_register(app, "slosh",           "fclaw_options.ini");
    sloshclawpatch_opt =    fclaw_clawpatch_2d_options_register(app, "slosh-clawpatch", "fclaw_options.ini");
    slosh_geo_opt =              fc2d_geoclaw_options_register(app, "slosh-geoclaw",   "fclaw_options.ini");
    slosh_user_opt =                    slosh_options_register(app, "slosh-user",      "fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    int first_arg;
    fclaw_exit_type_t vexit;
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);
    
        fclaw_global_t *radial_glob = fclaw_global_new_comm(mpicomm, size, rank);

        fclaw_options_store           (radial_glob, radial_gparms);
        fclaw_clawpatch_options_store (radial_glob, radial_clawpatchopt);
        fc2d_geoclaw_options_store      (radial_glob, radial_geoclawopt);
        radial_options_store            (radial_glob, radial_user_opt);

        radial_create_domain(radial_glob);

        /* Run the program */
        radial_run_program(radial_glob);

        fclaw_global_destroy(radial_glob);

    
        /* Create global structure which stores the domain, timers, etc */
        fclaw_global_t *slosh_glob = fclaw_global_new_comm(mpicomm, size, rank);

        /* Store option packages in glob */
        fclaw_options_store           (slosh_glob, slosh_fclaw_opt);
        fclaw_clawpatch_options_store (slosh_glob, sloshclawpatch_opt);
        fc2d_geoclaw_options_store      (slosh_glob, slosh_geo_opt);
        slosh_options_store             (slosh_glob, slosh_user_opt);

        slosh_create_domain(slosh_glob);

        slosh_run_program(slosh_glob);
        
        fclaw_global_destroy(slosh_glob);
    
    }

    fclaw_app_destroy (app);

    return 0;
}
