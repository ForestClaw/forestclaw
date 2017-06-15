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

#include "bowl_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include <fc2d_geoclaw.h>
#include <fc2d_geoclaw_options.h>

static int s_user_options_package_id = -1;

static void *
bowl_register (user_options_t* user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 0 = nomap; 1 = brick [0]");

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
bowl_check (user_options_t *user)
{
    if (user->example < 0 || user->example > 1)
    {
        fclaw_global_essentialf ("Option --user:example must be 0 or 1\n");
        return FCLAW_EXIT_QUIET;
    }
    return FCLAW_NOEXIT;
}

static void
bowl_destroy(user_options_t *user)
{
    /* Nothing to destroy */
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (user_options_t*) package;

    return bowl_register(user,opt);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return bowl_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (user_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    bowl_destroy (user);

    FCLAW_FREE (user);
}



static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    options_check,
    options_destroy
};

/* ------------- User options access functions --------------------- */

static
user_options_t* bowl_options_register (fclaw_app_t * app,
                                       const char *configfile)
{
    user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (user_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}


user_options_t* bowl_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}


static 
void bowl_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}


static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, 
                                fclaw_options_t* gparms, 
                                user_options_t* user)
{

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    /* Map unit square to disk using mapc2m_disk.f */
    int mi,mj;

    mi = gparms->mi;
    mj = gparms->mj;
    int a = 0; /* non-periodic */
    int b = 0;

    switch (user->example) {
    case 0:
        /* Size is set by [ax,bx] x [ay, by], set in .ini file */
        conn = p4est_connectivity_new_unitsquare();
        cont = fclaw2d_map_new_nomap();
        break;

    case 1:
        /* Square brick domain */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick(conn,mi,mj);
        cont = fclaw2d_map_new_nomap_brick(brick);
        break;

    default:
        SC_ABORT_NOT_REACHED ();
    }

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    return domain;
}

static
void run_program(fclaw2d_global_t* glob)
{
    fclaw2d_domain_t    **domain = &glob->domain;

    fclaw2d_domain_data_new(*domain);

    fclaw2d_vtable_initialize();
    fclaw2d_diagnostics_vtable_initialize();

    fc2d_geoclaw_vtable_initialize();

    bowl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fc2d_geoclaw_setup(glob);

    fclaw2d_initialize(glob);
    fclaw2d_run(glob);
    
    fc2d_geoclaw_finalize(glob);
    fclaw2d_finalize(glob);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                *options;
    fclaw_options_t             *gparms;
    fclaw2d_clawpatch_options_t *clawpatchopt;
    fc2d_geoclaw_options_t      *geoclawopt;
    user_options_t              *user;

    sc_MPI_Comm mpicomm;
    fclaw2d_domain_t* domain;
    fclaw2d_global_t* glob;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    gparms                   = fclaw_options_register(app,"fclaw_options.ini");
    clawpatchopt = fclaw2d_clawpatch_options_register(app, "fclaw_options.ini");
    geoclawopt        = fc2d_geoclaw_options_register(app, "fclaw_options.ini");
    user                      =  bowl_options_register(app,"fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!retval & !vexit)
    {
        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, gparms, user);
    
        glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(glob, domain);

        fclaw2d_options_store           (glob, gparms);
        fclaw2d_clawpatch_options_store (glob, clawpatchopt);
        fc2d_geoclaw_options_store      (glob, geoclawopt);
        bowl_options_store              (glob, user);

        /* Run the program */
        run_program(glob);

        fclaw2d_global_destroy(glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
