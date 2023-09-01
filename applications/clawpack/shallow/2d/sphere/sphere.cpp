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

#include "sphere_user.h"

#include <fclaw_forestclaw.h>
#include <fclaw_clawpatch.hpp>

#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>
#include <p4est_connectivity.h>

#include <fc2d_clawpack46.h>


typedef struct user_options
{
    int example;

    const char* latitude_string;
    double *latitude;

    const char* longitude_string;
    double *longitude;

    int is_registered;

} user_options_t;

static void *
options_register_user (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t* user = (user_options_t*) package;

    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 0,1 = cubedsphere; 2 = latlong [0]");

    fclaw_options_add_double_array(opt, 0, "latitude", &user->latitude_string,
                                   "-50 50", &user->latitude, 2,
                                   "[user] Latitude range (degrees) [-50 50]");

    fclaw_options_add_double_array(opt, 0, "longitude", &user->longitude_string,
                                   "0 360", &user->longitude, 2,
                                   "[user] Longitude range (degrees) [0 360]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
options_postprocess_user (fclaw_app_t * a, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;

    if (user->example == 0)
    {
        fclaw_options_convert_double_array (user->latitude_string, &user->latitude,2);
        fclaw_options_convert_double_array (user->longitude_string, &user->longitude,2);
    }
    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
options_check_user (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;
    if (user->example < 0 || user->example > 2) {
        fclaw_global_essentialf ("Option --user:example must be 0, 1, 2, 3 or 4\n");
        return FCLAW_EXIT_QUIET;
    }
    return FCLAW_NOEXIT;
}

static void
options_destroy_user (fclaw_app_t * a, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;
    /* Destroy arrays used in options  */
    if (user->example == 0)
    {
        fclaw_options_destroy_array((void*) user->latitude);
        fclaw_options_destroy_array((void*) user->longitude);
    }
}


static const
fclaw_app_options_vtable_t options_vtable_user =
{
    options_register_user,
    options_postprocess_user,
    options_check_user,
    options_destroy_user
};

static
void register_user_options (fclaw_app_t * app,
                            const char *configfile,
                            user_options_t* user)
{
    FCLAW_ASSERT (app != NULL);

    /* sneaking the version string into the package pointer */
    /* when there are more parameters to pass, create a structure to pass */
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);
}

static
void run_program(fclaw_app_t* app)
{
    sc_MPI_Comm            mpicomm;

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t	     *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    fclaw_options_t   *gparms;
    user_options_t  *user;

    /* Used locally */
    int mi, mj, a,b;
    double pi = M_PI;
    double rotate[2];

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

    gparms = fclaw_forestclaw_get_options(app);
    user = (user_options_t*) fclaw_app_get_user(app);

    /* ---------------------------------------------------------------
       Mapping geometry
       --------------------------------------------------------------- */
    mi = gparms->mi;
    mj = gparms->mj;
    a = gparms->periodic_x;
    b = gparms->periodic_y;
    rotate[0] = pi*gparms->theta/180.0;
    rotate[1] = pi*gparms->phi/180.0;

    switch (user->example)
    {
    case 0:
    case 1:
        /* Cubed sphere */
        conn = p4est_connectivity_new_cubed();
        cont = fclaw2d_map_new_cubedsphere(gparms->scale,gparms->shift,rotate);
        break;
    case 2:
        /* latlong */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        cont = fclaw2d_map_new_latlong(brick,gparms->scale,
                                       user->latitude,user->longitude,a,b);
        break;
    default:
        SC_ABORT_NOT_REACHED ();
    }

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_set_app (domain,app);

    fclaw_domain_data_new(domain);

    sphere_link_solvers(domain);

    fclaw_initialize(&domain);
    fclaw_run(&domain);
    fclaw_finalize(&domain);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    user_options_t    suser_options, *user = &suser_options;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, user);

    /* Register packages */
    fclaw_forestclaw_register(app,"fclaw_options_ini");
    fc2d_clawpack46_register(app,"fclaw_options.ini");

    /* User defined options (defined above) */
    register_user_options (app, "fclaw_options.ini", user);

    /* Read configuration file(s) */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Link packages to patches */

    fclaw2d_clawpatch_link_app(app);

    if (!vexit)
    {
        run_program(app);
    }

    fclaw_forestclaw_destroy(app);
    fclaw_app_destroy (app);

    return 0;
}
