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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>

#include "metric_user.h"


static void *
options_register_user (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t* user = (user_options_t*) package;

    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "1 cart; 2 5-patch square; 3 squared-disk; "    \
                        "4 pillow disk; 5 pillowdisk5; 6 pillow sphere " \
                        "7 cubed sphere; 8 torus");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 0.4,
                           "[user] Ratio used for squared- and pillow-disk [0.4]");

    sc_options_add_double (opt, 0, "beta", &user->beta, 0.4,
                           "[user] Ratio of inner to outer radii in torus [0.4]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 4,
                           "Clawpack_version (4 or 5) [4]");

    user->is_registered = 1;
    return NULL;
}


static fclaw_exit_type_t
options_check_user (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;

    if (user->example < 0 || user->example > 8) {
        fclaw_global_essentialf ("Option --user:example must be 1-8\n");
        return FCLAW_EXIT_QUIET;
    }
    return FCLAW_NOEXIT;
}

static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register_user,
    NULL,
    options_check_user,
    NULL
};

static
void register_user_options (fclaw_app_t * app,
                            const char *configfile,
                            user_options_t* user)
{
    FCLAW_ASSERT (app != NULL);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);
}

const user_options_t* metric_user_get_options(fclaw2d_domain_t* domain)
{
    fclaw_app_t* app;
    app = fclaw2d_domain_get_app(domain);

    const user_options_t* user = (user_options_t*) fclaw_app_get_user(app);

    return (user_options_t*) user;
}

void run_program(fclaw_app_t* app)
{
    sc_MPI_Comm            mpicomm;

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t	     *domain;
    fclaw_map_context_t    *cont = NULL;

    fclaw_options_t* gparms;
    user_options_t suser_options, *user = &suser_options;

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

    gparms = fclaw_forestclaw_get_options(app);
    user = (user_options_t*) fclaw_app_get_user(app);

    /* ---------------------------------------------------------------
       Domain geometry
       -------------------------------------------------------------- */
    double pi = M_PI;
    double rotate[2];

    rotate[0] = pi*gparms->phi/180.0;
    rotate[1] = pi*gparms->theta/180.0;

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

    switch (user->example) {
    case 0:
    case 1:
        /* Map [0,1]x[0,1] to [-1,1],[-1,1] */
        conn = p4est_connectivity_new_unitsquare();
        cont = fclaw2d_map_new_cart (gparms->scale, gparms->shift,rotate);
        break;
    case 2:
        /* Map [0,1],[0,1] to five patch square [-1,1]x[-1,1] */
        conn = p4est_connectivity_new_disk ();
        cont = fclaw2d_map_new_fivepatch (gparms->scale,gparms->shift,
                                          rotate, user->alpha);
        break;
    case 3:
        /* Map [0,1],[0,1] to squared disk (radius 1, centered at the origin) */
        conn = p4est_connectivity_new_disk ();
        cont = fclaw2d_map_new_squareddisk (gparms->scale,gparms->shift,
                                            rotate,user->alpha);
        break;
    case 4:
        /* Map [0,1],[0,1] to pillow disk */
        conn = p4est_connectivity_new_unitsquare ();
        cont = fclaw2d_map_new_pillowdisk (gparms->scale,gparms->shift,rotate);
        break;
    case 5:
        /* Map [0,1]x[0,1] to five patch --> pillow disk */
        conn = p4est_connectivity_new_disk ();
        cont = fclaw2d_map_new_pillowdisk5 (gparms->scale,gparms->shift,
                                            rotate,user->alpha);
        break;
    case 6:
        /* Map [0,1]x[0,1] to five patch --> pillow disk */
        conn = p4est_connectivity_new_pillow ();
        cont = fclaw2d_map_new_pillowsphere (gparms->scale,gparms->shift,rotate);
        break;
    case 7:
      /* Map [0,1]x[0,1] to five patch --> pillow disk */
        conn = p4est_connectivity_new_cubed ();
        cont = fclaw2d_map_new_cubedsphere (gparms->scale,gparms->shift,rotate);
        break;
    case 8:
        /* Map [0,1]x[0,1] to five patch --> pillow disk */
        conn = p4est_connectivity_new_periodic ();
        cont = fclaw2d_map_new_torus (gparms->scale,gparms->shift,
                                      rotate,user->beta);
        break;
    default:
        SC_ABORT_NOT_REACHED ();
  }

  domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

  /* ---------------------------------------------------------- */

  fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_INFO);
  fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);


  /* ---------------------------------------------------------------
     Set domain data.
     --------------------------------------------------------------- */
  fclaw2d_domain_data_new(domain);
  fclaw2d_domain_set_app(domain,app);

  /* Link other routines that need to be included. */
  metric_link_patch(domain);

  /* --------------------------------------------------
     Initialize and run the simulation
     -------------------------------------------------- */

  fclaw2d_initialize(&domain);
  fclaw2d_run(&domain);
  fclaw2d_finalize(&domain);
}

int
main (int argc, char **argv)
{
  fclaw_app_t *app;
  int first_arg;
  fclaw_exit_type_t vexit;

  /* Options */
  user_options_t suser_options, *user = &suser_options;

  /* Initialize application */
  app = fclaw_app_new (&argc, &argv, user);

  fclaw_forestclaw_register(app,"fclaw_options.ini");

  /* User defined options (defined above) */
  register_user_options (app, "fclaw_options.ini", user);
  fc2d_clawpack46_register(app,"fclaw_options.ini");    /* [clawpack46] */


  /* Read configuration file(s) */
  vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

  /* No packages to register */
  fclaw2d_clawpatch_link_app(app);


  if (!vexit)
  {
      run_program(app);
  }

  fclaw_forestclaw_destroy(app);
  fclaw_app_destroy (app);

  return 0;
}
