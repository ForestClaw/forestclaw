/*
  Copyright (c) 2019-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "heat/heat_user.h"
#include "phasefield/phasefield_user.h"
    
#include <fclaw2d_include_all.h>

#include <fclaw2d_output.h>
#include <fclaw2d_diagnostics.h>

#include <fclaw2d_elliptic_solver.h>

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.h>

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_options.h>


int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    fclaw_options_t             *heat_fclaw_opt;

    fclaw_clawpatch_options_t *heat_clawpatch_opt;
    fc2d_thunderegg_options_t   *heat_mg_opt;
    heat_options_t              *heat_user_opt;

    fclaw_global_t            *heat_glob;
    fclaw_domain_t            *heat_domain;

    fclaw_options_t             *phasefield_fclaw_opt;

    fclaw_clawpatch_options_t *phasefield_clawpatch_opt;
    fc2d_thunderegg_options_t   *phasefield_mg_opt;
    phasefield_options_t        *phasefield_user_opt;

    fclaw_global_t            *phasefield_glob;
    fclaw_domain_t            *phasefield_domain;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    heat_fclaw_opt =                   fclaw_options_register(app, "heat",            "fclaw_options.ini");
    heat_clawpatch_opt =   fclaw_clawpatch_options_register_2d(app, "heat-clawpatch",  "fclaw_options.ini");
    heat_mg_opt =            fc2d_thunderegg_options_register(app, "heat-thunderegg", "fclaw_options.ini");
    heat_user_opt =                     heat_options_register(app, "heat-user",       "fclaw_options.ini");  

    phasefield_fclaw_opt =                   fclaw_options_register(app, "phasefield",           "fclaw_options.ini");
    phasefield_clawpatch_opt =   fclaw_clawpatch_options_register_2d(app, "phasefield-clawpatch",  "fclaw_options.ini");
    phasefield_mg_opt =            fc2d_thunderegg_options_register(app, "phasefield-thunderegg", "fclaw_options.ini");
    phasefield_user_opt =               phasefield_options_register(app, "phasefield-user",       "fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        heat_domain = heat_create_domain(mpicomm, heat_fclaw_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        heat_glob = fclaw_global_new();
        fclaw_global_store_domain(heat_glob, heat_domain);

        /* Store option packages in glob */
        fclaw2d_options_store           (heat_glob, heat_fclaw_opt);
        fclaw_clawpatch_options_store (heat_glob, heat_clawpatch_opt);
        fc2d_thunderegg_options_store    (heat_glob, heat_mg_opt);
        heat_options_store            (heat_glob, heat_user_opt);

        heat_run_program(heat_glob);

        fclaw_global_destroy(heat_glob);        
        phasefield_domain = phasefield_create_domain(mpicomm, phasefield_fclaw_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        phasefield_glob = fclaw_global_new();
        fclaw_global_store_domain(phasefield_glob, phasefield_domain);

        /* Store option packages in glob */
        fclaw2d_options_store           (phasefield_glob, phasefield_fclaw_opt);
        fclaw_clawpatch_options_store (phasefield_glob, phasefield_clawpatch_opt);
        fc2d_thunderegg_options_store    (phasefield_glob, phasefield_mg_opt);
        phasefield_options_store            (phasefield_glob, phasefield_user_opt);

        phasefield_run_program(phasefield_glob);

        fclaw_global_destroy(phasefield_glob); 
    }
    
    fclaw_app_destroy (app);

    return 0;
}
