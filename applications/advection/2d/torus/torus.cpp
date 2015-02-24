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

#include <forestclaw2d.h>
#include <fclaw_options.h>

#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw_options.h>

#include <p4est_connectivity.h>

#include <fc2d_clawpack46.H>

#include "torus_user.H"

typedef struct user_options
{
    int example;
    double alpha;
    double beta;

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
                        "[user] 1 = cart; 2 = torus; 3 = lat-long; 4 = annulus [2]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 0.4,
                           "[user] Ratio r/R, r=outer radius, R=inner radius " \
                           "(used for torus) [0.4]");

    fclaw_options_add_double_array(opt, 0, "latitude", &user->latitude_string,
                                   "-50 50", &user->latitude, 2,
                                   "[user] Latitude range (degrees) [-50 50]");

    fclaw_options_add_double_array(opt, 0, "longitude", &user->longitude_string,
                                   "0 360", &user->longitude, 2,
                                   "[user] Longitude range (degrees) [0 360]");

    sc_options_add_double (opt, 0, "beta", &user->beta, 0.4,
                           "[user] Inner radius of annulus [0.4]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
options_postprocess_user (fclaw_app_t * a, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;

    if (user->example == 3)
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
    if (user->example < 0 || user->example > 4) {
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
    if (user->example == 3)
    {
        fclaw_options_destroy_array((void*) user->latitude);
        fclaw_options_destroy_array((void*) user->longitude);
    }
}


static const fclaw_app_options_vtable_t options_vtable_user =
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

void run_program(fclaw_app_t* app, amr_options_t* gparms,
                 fc2d_clawpack46_options_t* clawpack_options,
                 user_options_t* user)
{
    sc_MPI_Comm            mpicomm;

    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t	     *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    /* Used locally */
    double pi = M_PI;
    double rotate[2];
    int mi, mj, a,b;

    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

    /* ---------------------------------------------------------------
       Mapping geometry
       --------------------------------------------------------------- */
    mi = gparms->mi;
    mj = gparms->mj;
    rotate[0] = pi*gparms->theta/180.0;
    rotate[1] = pi*gparms->phi/180.0;
    a = gparms->periodic_x;
    b = gparms->periodic_y;

    switch (user->example)
    {
    case 0:
        FCLAW_ASSERT(mi == 1 && mj == 1);  /* assumes square domain */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        cont = fclaw2d_map_new_nomap();
        break;
    case 1:
        /* Cartesian [-1,1]x[-1,1] */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick(conn,mi,mj);
        cont = fclaw2d_map_new_cart(brick, gparms->scale, gparms->shift, rotate);
        break;
    case 2:
        /* torus */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick(conn,mi,mj);
        cont = fclaw2d_map_new_torus(brick,gparms->scale,gparms->shift,rotate,user->alpha);
        break;
    case 3:
        /* Lat-long example */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick(conn,mi,mj);
        cont = fclaw2d_map_new_latlong(brick,gparms->scale,gparms->shift,
                                       rotate,user->latitude,user->longitude,a,b);
        break;
    case 4:
        /* Annulus */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick(conn,mi,mj);
        cont = fclaw2d_map_new_annulus(brick,gparms->scale,gparms->shift,
                                       rotate,user->beta);
        break;

    default:
        SC_ABORT_NOT_REACHED ();
    }

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    init_domain_data(domain);

    set_domain_parms(domain,gparms);
    fc2d_clawpack46_set_options(domain,clawpack_options);

    link_problem_setup(domain,fc2d_clawpack46_setprob);

    torus_link_solvers(domain);

    link_regrid_functions(domain,
                          torus_patch_tag4refinement,
                          torus_patch_tag4coarsening);
    amrinit(&domain);
    amrrun(&domain);
    amrreset(&domain);

    fclaw2d_map_destroy(cont);
}


int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t              *options;
    amr_options_t             samr_options, *gparms = &samr_options;
    fc2d_clawpack46_options_t  sclawpack_options, *clawpack_options = &sclawpack_options;
    user_options_t                suser_options, *user = &suser_options;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, user);
    options = fclaw_app_get_options (app);
    fclaw_package_container_new(app);

    fclaw_options_register_general (app, "fclaw_options.ini", gparms);
    fc2d_clawpack46_options_register(app,"fclaw_options_ini",clawpack_options);

    /* User defined options (defined above) */
    register_user_options (app, "fclaw_options.ini", user);


    /* Read configuration file(s) */
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Link packages to patches */
    fc2d_clawpack46_package_register(app,clawpack_options);
    link_app_to_clawpatch(app);

    if (!retval & !vexit)
    {
        run_program(app, gparms, clawpack_options,user);
    }

    fclaw_package_container_destroy(app);
    fclaw_app_destroy (app);

    return 0;
}
