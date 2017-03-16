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

#include "torus_user.h"

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>

#include <fclaw2d_map.h>
#include <fclaw2d_map_brick.h>
#include <fclaw2d_map_query.h>
#include <p4est_connectivity.h>

#include <fc2d_clawpack46.h>

static int s_user_package_id = -1;

static void *
options_register_user (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t* user = (user_options_t*) package;

    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 0 = torus; 1 = twisted torus [0]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                        "[user] Clawpack version (4 or 5) [5]");


    sc_options_add_double (opt, 0, "alpha", &user->alpha, 0.4,
                           "[user] Ratio r/R, r=outer radius, R=inner radius " \
                           "(used for torus) [0.4]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
options_check_user (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t* user = (user_options_t*) package;

    if (user->example < 0 || user->example > 1)
    {
        fclaw_global_essentialf
            ("Option --user:example must be 0 or 1\n");
        return FCLAW_EXIT_QUIET;
    }
    return FCLAW_NOEXIT;

}

static const
fclaw_app_options_vtable_t options_vtable_user =
{
    options_register_user,
    NULL,
    options_check_user,
    NULL,
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
void user_set_options (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_package_id = id;
}

const user_options_t* torus_user_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_package_id;
    return (user_options_t*) 
            fclaw_package_get_options(glob, id);    
}

static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, amr_options_t* gparms, user_options_t* user)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    /* Used locally */
    double pi = M_PI;
    double rotate[2];
    int mi, mj, a,b;    

    /* ---------------------------------------------------------------
       Mapping geometry
       --------------------------------------------------------------- */
    mi = gparms->mi;
    mj = gparms->mj;
    rotate[0] = pi*gparms->theta/180.0;
    rotate[1] = 0;  /* Don't rotate through phi */

    a = 1;  /* Torus is periodic in both directions */
    b = 1;

    /* torus */
    conn  = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick(conn,mi,mj);
    cont  = fclaw2d_map_new_torus(brick,gparms->scale,gparms->shift,rotate,
                                  user->alpha,user->example);

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);

    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
    return domain;
}
static
void run_program(fclaw2d_global_t* glob)
{
    user_options_t  *user;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    user = (user_options_t*) torus_user_get_options(glob);

    if (user->claw_version == 4)
    {
      fc2d_clawpack46_set_vtable_defaults();
    }
    else if (user->claw_version == 5)
    {
      fc2d_clawpack5_set_vtable_defaults();
    }

    torus_link_solvers(glob);

    fclaw2d_initialize(glob);
    fclaw2d_run(glob);
    fclaw2d_finalize(glob);
}


int
main (int argc, char **argv)
{
    int first_arg;
    fclaw_app_t *app;
    fclaw2d_global_t *glob;
    fclaw_exit_type_t vexit;

    sc_MPI_Comm mpicomm;
    fclaw2d_domain_t *domain;

    /* Options */
    sc_options_t                *options;
    user_options_t              suser, *user = &suser;
    amr_options_t               *gparms;
    fclaw2d_clawpatch_options_t *clawpatchopt;
    fc2d_clawpack46_options_t   *claw46opt;
    fc2d_clawpack5_options_t    *claw5opt;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, user);

    /* Register packages */
    gparms = fclaw2d_forestclaw_options_register(app,"fclaw_options.ini");
    clawpatchopt = fclaw2d_clawpatch_options_register(app, "fclaw_options.ini");
    claw46opt = fc2d_clawpack46_options_register(app,"fclaw_options.ini");
    claw5opt = fc2d_clawpack5_options_register(app,"fclaw_options.ini");
    register_user_options(app,"fclaw_options.ini",user);  /* [user] */
    /* User defined options (defined above) */
    register_user_options (app, "fclaw_options.ini", user);

    /* Read configuration file(s) */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* at this point gparms is valid */
    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
    domain = create_domain(mpicomm, gparms, user);

    glob = fclaw2d_global_new();
    fclaw2d_global_set_domain(glob, domain);
    
    fclaw2d_forestclaw_set_options (glob, gparms);
    fclaw2d_clawpatch_set_options (glob, clawpatchopt);
    fc2d_clawpack46_set_options (glob, claw46opt);
    fc2d_clawpack5_set_options (glob, claw5opt);
    user_set_options (glob, user);


    if (!retval & !vexit)
    {
        run_program(glob);
    }

    fclaw2d_global_destroy (glob);
    fclaw_app_destroy (app);

    return 0;
}
