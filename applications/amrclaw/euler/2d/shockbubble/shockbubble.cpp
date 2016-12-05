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

#include "shockbubble_user.h"

#include "fc2d_clawpack46.h"
#include "fc2d_clawpack5.h"

static void *
options_register_user (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t* user = (user_options_t*) package;

    /* [user] User options */
    sc_options_add_int (opt, 0, "idisc", &user->idisc,    2, "[user] idisc [2]");
    sc_options_add_double (opt, 0, "gamma", &user->gamma, 1.4, "[user] gamma [1.4]");

    sc_options_add_double (opt, 0, "x0",    &user->x0,    0.5, "[user] x0 [0.5]");
    sc_options_add_double (opt, 0, "y0",    &user->y0,    0.0, "[user] y0 [0.0]");
    sc_options_add_double (opt, 0, "r0",    &user->r0,    0.2, "[user] r0 [0.2]");
    sc_options_add_double (opt, 0, "rhoin", &user->rhoin, 0.1, "[user] rhoin [0.1]");
    sc_options_add_double (opt, 0, "pinf",  &user->pinf,  5.0, "[user] pinf [5.0x]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                           "Clawpack_version (4 or 5) [5]");

    user->is_registered = 1;
    return NULL;
}

static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register_user,
    NULL,
    NULL,
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

user_options_t* shockbubble_user_get_options(fclaw2d_domain_t* domain)
{
    const user_options_t* user = (user_options_t*) fclaw2d_domain_get_user_options(domain);

    return (user_options_t*) user;
}


void run_program(fclaw_app_t* app)
{
    sc_MPI_Comm            mpicomm;

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t	     *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    amr_options_t *gparms;

    int mi, mj, a,b;

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
    gparms = fclaw_forestclaw_get_options(app);

    mi = gparms->mi;
    mj = gparms->mj;
    a = gparms->periodic_x;
    b = gparms->periodic_y;

    /* Use [ax,bx]x[ay,by] */

    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick(conn,mi,mj);
    cont = fclaw2d_map_new_nomap_brick(brick);

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

    /* ---------------------------------------------------------- */
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_INFO);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    fclaw2d_domain_data_new(domain);
    fclaw2d_domain_set_app(domain,app);

    shockbubble_link_solvers(domain);

    fclaw2d_initialize(&domain);
    fclaw2d_run(&domain);
    fclaw2d_finalize(&domain);

    fclaw2d_map_destroy(cont);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                  *options;
    user_options_t suser, *user = &suser;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, user);

    /* Register packages */
    fclaw_forestclaw_register(app,"fclaw_options.ini");
    fc2d_clawpack46_register(app,"fclaw_options.ini");
    fc2d_clawpack5_register(app,"fclaw_options.ini");

    register_user_options(app,"fclaw_options.ini",user);

    /* Read configuration file(s) */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    fclaw2d_clawpatch_link_app(app);

    if (!retval & !vexit)
    {
        run_program(app);
    }

    fclaw_forestclaw_destroy(app);
    fclaw_app_destroy (app);

    return 0;
}
