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

#include <amr_single_step.h>
#include <fclaw2d_map.h>
#include <p4est_connectivity.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw_options.h>

#include <fclaw2d_map_query.h>

#include <fclaw_package.h>
#include <fclaw_register.h>
#include <fc2d_clawpack46.H>

#include "sphere_user.H"

typedef struct user_options
{
    int example;
    int is_registered;

} user_options_t;

static void *
    options_register_user (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t* user = (user_options_t*) package;

    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 1 for pillow grid, "    \
                        "2 for cubed sphere ");
    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
    options_check_user (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;
    if (user->example < 1 || user->example > 2) {
        fclaw_global_essentialf ("Option --user:example must be 1 or 2\n");
        return FCLAW_EXIT_ERROR;
    }
    return FCLAW_NOEXIT;
}

static void
    options_destroy_user(fclaw_app_t* app, void* package, void *registered)
{
    user_options_t* user = (user_options_t*) package;
    FCLAW_ASSERT(user != NULL);
    FCLAW_FREE(user);
};

static const fclaw_app_options_vtable_t options_vtable_user = {
    options_register_user,
    NULL,      /* options_postprocess_user */
    options_check_user,
    options_destroy_user       /* options_destroy_user */
};


void static
    register_user_options (fclaw_app_t * app,
                           const char *configfile)
{
    user_options_t* user;

    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC(user_options_t,1);

    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);
    fclaw_app_set_attribute(app,"user",user);
}

static
    void run_program(fclaw_app_t* app)
    {
        sc_MPI_Comm            mpicomm;

        /* Mapped, multi-block domain */
        p4est_connectivity_t     *conn = NULL;
        fclaw2d_domain_t	 *domain;
        fclaw2d_map_context_t    *cont = NULL;

        amr_options_t              *gparms;
        fc2d_clawpack46_options_t  *clawpack_options;
        user_options_t             *user;

        /* Used locally */
        double pi = M_PI;
        double rotate[2];

        gparms = fclaw_forestclaw_get_options(app);
        clawpack_options = fc2d_clawpack46_get_options(app);
        user = (user_options_t*) fclaw_app_get_attribute(app,"user",NULL);

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

        rotate[0] = pi*gparms->theta/180.0;
        rotate[1] = pi*gparms->phi/180.0;

        switch (user->example) {
        case 1:
            conn = p4est_connectivity_new_pillow();
            cont = fclaw2d_map_new_pillowsphere(gparms->scale,gparms->shift,rotate);
            break;
        case 2:
            conn = p4est_connectivity_new_cubed();
            cont = fclaw2d_map_new_cubedsphere(gparms->scale,gparms->shift,rotate);
            break;
        default:
            SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
        }

        domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

        fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
        fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

        /* ---------------------------------------------------------------
           Set domain data.
           --------------------------------------------------------------- */
        init_domain_data(domain);

        /* Store parameters */
        set_domain_parms(domain,gparms);
        fc2d_clawpack46_set_options (domain,clawpack_options);

        /* Link solvers to the domain */
        link_problem_setup(domain,fc2d_clawpack46_setprob);

        sphere_link_solvers(domain);

        /* --------------------------------------------------
           Initialize and run the simulation
           -------------------------------------------------- */
        amrinit(&domain);
        amrrun(&domain);
        amrreset(&domain);

        /* --------------------------------------------------
           Clean up the mapping context.
           -------------------------------------------------- */
        fclaw2d_map_destroy (cont);
    }


int main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;
    sc_options_t             *options;
    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    fclaw_forestclaw_register(app,"fclaw_options.ini");
    fc2d_clawpack46_register(app,"fclaw_options.ini");

    /* User options */
    register_user_options (app, "fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    link_app_to_clawpatch(app);

    /* Run program */
    if (!retval & !vexit)
    {
        run_program(app);
    }

    fclaw_forestclaw_destroy(app);
    fclaw_app_destroy (app);

    return 0;
}
