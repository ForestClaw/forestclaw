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

#include "filament_user.h"

#include "../all/advection_user.h"

static void
store_domain_map (fclaw2d_global_t * glob, fclaw_options_t * fclaw_opt,
                  user_options_t * user,
                  fclaw2d_clawpatch_options_t * clawpatch_opt)
{
    /* Mapped, multi-block domain */
    fclaw2d_domain_t *domain = NULL;
    fclaw2d_map_context_t *cont = NULL, *brick = NULL;

    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;
    int a = 0; /* non-periodic */
    int b = 0;

    int mx = clawpatch_opt->mx;
    int minlevel = fclaw_opt->minlevel;
    
    switch (user->example) {
    case 1:
        /* Square brick domain */
        domain =
            fclaw2d_domain_new_brick (glob->mpicomm, mi, mj, a, b,
                                      fclaw_opt->minlevel);
        brick = fclaw2d_map_new_brick (domain, mi, mj, a, b);
        
        /* Square in [-1,1]x[-1,1], shifted by (1,1,0) */
        cont = fclaw2d_map_new_cart(brick,
                                    fclaw_opt->scale,
                                    fclaw_opt->shift);
        break;
    case 2:
        if (mi*mx*pow_int(2,minlevel) < 32)
        {
            fclaw_global_essentialf("The five patch mapping requires mi*mx*2^minlevel > 32\n");
            exit(0);

        }
        /* Five patch square domain */
        domain =
            fclaw2d_domain_new_disk (glob->mpicomm, 0, 0,
                                     fclaw_opt->minlevel);
        cont =
            fclaw2d_map_new_fivepatch (fclaw_opt->scale, fclaw_opt->shift,
                                       user->alpha);
        break;
    case 3:
        /* bilinear square domain : maps to [-1,1]x[-1,1] */
        FCLAW_ASSERT (mi == 2 && mj == 2);
        domain =
            fclaw2d_domain_new_brick (glob->mpicomm, mi, mj, a, b,
                                      fclaw_opt->minlevel);
        brick = fclaw2d_map_new_brick (domain, mi, mj, a, b);
        cont = fclaw2d_map_new_bilinear (brick, 
                                         fclaw_opt->scale,
                                         fclaw_opt->shift, 
                                         user->center);
        break;

    default:
        SC_ABORT_NOT_REACHED ();
    }

    fclaw2d_domain_list_levels (domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors (domain, FCLAW_VERBOSITY_DEBUG);
    fclaw2d_global_store_domain (glob, domain);
    fclaw2d_global_store_map (glob, cont);
}

static
void run_program(fclaw2d_global_t* glob)
{
    user_options_t             *user;

    user = (user_options_t*) filament_get_options(glob);


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
    fclaw2d_run(glob);
    fclaw2d_finalize(glob);
}



int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    user_options_t              *user_opt;
    fclaw_options_t             *fclaw_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;
    fc2d_clawpack46_options_t   *claw46_opt;
    fc2d_clawpack5_options_t    *claw5_opt;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    fclaw_opt                    = fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt    = fclaw2d_clawpatch_options_register(app, "clawpatch",  "fclaw_options.ini");
    claw46_opt         = fc2d_clawpack46_options_register(app, "clawpack46", "fclaw_options.ini");
    claw5_opt           = fc2d_clawpack5_options_register(app, "clawpack5",  "fclaw_options.ini");
    user_opt =                  filament_options_register(app,               "fclaw_options.ini");  

    /* Read configuration file(s) */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        /* Options have been checked and are valid */
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);
        fclaw2d_global_t *glob =
            fclaw2d_global_new_comm (mpicomm, size, rank);
        store_domain_map (glob, fclaw_opt, user_opt, clawpatch_opt);

        fclaw2d_options_store            (glob, fclaw_opt);
        fclaw2d_clawpatch_options_store  (glob, clawpatch_opt);
        fc2d_clawpack46_options_store    (glob, claw46_opt);
        fc2d_clawpack5_options_store     (glob, claw5_opt);
        filament_options_store           (glob, user_opt);

        run_program(glob);

        fclaw2d_global_destroy(glob);        
    }

    fclaw_app_destroy (app);

    return 0;
}
