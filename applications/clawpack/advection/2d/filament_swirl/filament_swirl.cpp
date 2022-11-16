/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "filament/filament_user.h"
#include "swirl/swirl_user.h"
#include "user.h"

#include "../all/advection_user.h"

#include <fclaw_filesystem.h>
#include <fclaw_base.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_global.h>

static
void filament_initialize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    filament_options_t             *user;

    user = (filament_options_t*) filament_get_options(glob);


    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    if (user->claw_version == 4)
    {
      fc2d_clawpack46_solver_initialize(glob);
    }
    else if (user->claw_version == 5)
    {
      fc2d_clawpack5_solver_initialize(glob);
    }

    filament_link_solvers(glob);

    fclaw2d_initialize(glob);

    fclaw2d_clear_global_context(glob);
}

static
void filament_finalize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    fclaw2d_problem_setup(glob);
    fclaw2d_finalize(glob);

    fclaw2d_clear_global_context(glob);
}
static
void swirl_initialize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    const swirl_options_t           *swirl_opt;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    swirl_opt = swirl_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (swirl_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize(glob);
    }
    else if (swirl_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize(glob);
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw2d_initialize(glob);

    fclaw2d_clear_global_context(glob);
}
static
void swirl_finalize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    fclaw2d_problem_setup(glob);
    fclaw2d_finalize(glob);

    fclaw2d_clear_global_context(glob);
}

static
void run_programs(fclaw2d_global_t* globs[])
{
    user_run(globs,2);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                *options;

    filament_options_t          *filament_user_opt;
    fclaw_options_t             *filament_fclaw_opt;
    fclaw2d_clawpatch_options_t *filament_clawpatch_opt;
    fc2d_clawpack46_options_t   *filament_claw46_opt;
    fc2d_clawpack5_options_t    *filament_claw5_opt;

    swirl_options_t             *swirl_user_opt;
    fclaw_options_t             *swirl_fclaw_opt;
    fclaw2d_clawpatch_options_t *swirl_clawpatch_opt;
    fc2d_clawpack46_options_t   *swirl_claw46_opt;
    fc2d_clawpack5_options_t    *swirl_claw5_opt;

    fclaw2d_global_t         *filament_glob;
    fclaw2d_domain_t         *filament_domain;

    fclaw2d_global_t         *swirl_glob;
    fclaw2d_domain_t         *swirl_domain;

    sc_MPI_Comm mpicomm;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    filament_fclaw_opt                    = fclaw_options_register(app,  "filament",           "fclaw_options.ini");
    filament_clawpatch_opt    = fclaw2d_clawpatch_options_register(app, "filament-clawpatch",  "fclaw_options.ini");
    filament_claw46_opt         = fc2d_clawpack46_options_register(app, "filament-clawpack46", "fclaw_options.ini");
    filament_claw5_opt           = fc2d_clawpack5_options_register(app, "filament-clawpack5",  "fclaw_options.ini");
    filament_user_opt =                  filament_options_register(app, "filament-user",       "fclaw_options.ini");  

    swirl_fclaw_opt =                   fclaw_options_register(app, "swirl",            "fclaw_options.ini");
    swirl_clawpatch_opt =   fclaw2d_clawpatch_options_register(app, "swirl-clawpatch",  "fclaw_options.ini");
    swirl_claw46_opt =        fc2d_clawpack46_options_register(app, "swirl-clawpack46", "fclaw_options.ini");
    swirl_claw5_opt =          fc2d_clawpack5_options_register(app, "swirl-clawpack5",  "fclaw_options.ini");
    swirl_user_opt =                    swirl_options_register(app, "swirl-user",       "fclaw_options.ini");  

    /* Read configuration file(s) */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!retval & !vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

        /* Domains */
        filament_domain = filament_create_domain(mpicomm, filament_fclaw_opt, filament_user_opt,filament_clawpatch_opt);
        swirl_domain = swirl_create_domain(mpicomm, swirl_fclaw_opt);
            

        /* Globs */
        filament_glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(filament_glob, filament_domain);

        fclaw2d_options_store            (filament_glob, filament_fclaw_opt);
        fclaw2d_clawpatch_options_store  (filament_glob, filament_clawpatch_opt);
        fc2d_clawpack46_options_store    (filament_glob, filament_claw46_opt);
        fc2d_clawpack5_options_store     (filament_glob, filament_claw5_opt);
        filament_options_store           (filament_glob, filament_user_opt);

        swirl_glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(swirl_glob, swirl_domain);

        fclaw2d_options_store           (swirl_glob, swirl_fclaw_opt);
        fclaw2d_clawpatch_options_store (swirl_glob, swirl_clawpatch_opt);
        fc2d_clawpack46_options_store   (swirl_glob, swirl_claw46_opt);
        fc2d_clawpack5_options_store    (swirl_glob, swirl_claw5_opt);
        swirl_options_store             (swirl_glob, swirl_user_opt);

        /* initialize */
        filament_initialize(filament_glob);
        swirl_initialize(swirl_glob);

        /* run */
        fclaw2d_global_t* globs[2];
        globs[0] = filament_glob;
        globs[1] = swirl_glob;
        run_programs(globs);

        /* finalzie */
        filament_finalize(filament_glob);
        swirl_finalize(swirl_glob);


        /* destroy */
        fclaw2d_global_destroy(filament_glob);        
        fclaw2d_global_destroy(swirl_glob);

    }

    fclaw_app_destroy (app);

    return 0;
}
