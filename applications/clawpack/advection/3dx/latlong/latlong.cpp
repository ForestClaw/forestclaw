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

/* Needed for 2d mappings */
#include <fclaw2d_map.h>

#include "latlong_user.h"

#include "../all/advection_user.h"

static
fclaw_domain_t* create_domain(sc_MPI_Comm mpicomm, 
                                fclaw_options_t* fclaw_opt, 
                                user_options_t* user)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    /* ---------------------------------------------------------------
       Mapping geometry
       --------------------------------------------------------------- */
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;
    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    FCLAW_ASSERT(fclaw_opt->manifold != 0);

    /* Used locally */
    double pi = M_PI;
    double rotate[2];

    rotate[0] = pi*fclaw_opt->theta/180.0;
    rotate[1] = pi*fclaw_opt->phi/180.0;


    switch (user->example) 
    {
    case 1:
        /* Latlong */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        cont = fclaw2d_map_new_latlong(brick,fclaw_opt->scale,
                                       rotate,
                                       user->latitude, 
                                       user->longitude,
                                       a,b);
        break;
    default:
        SC_ABORT_NOT_REACHED ();
    }

    latlong_map_extrude(cont,user->maxelev);

    domain = fclaw_domain_wrap_2d(fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont));
    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
    return domain;
}

static
void run_program(fclaw_global_t* glob)
{
    const user_options_t           *user_opt;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    
    user_opt = latlong_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (user_opt->claw_version == 4)
    {
        fc3d_clawpack46_solver_initialize(glob);
    }
    else if (user_opt->claw_version == 5)
    {
        fclaw_global_essentialf("latlong : clawpack 5.0 not yet implemented with mappings\n");
        exit(0);
    }

    latlong_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw_initialize(glob);
    fclaw_run(glob);
    fclaw_finalize(glob);
}

int
main (int argc, char **argv)
{
    int first_arg;
    fclaw_app_t *app;
    fclaw_exit_type_t vexit;

    /* Options */
    user_options_t              *user_opt;
    fclaw_options_t             *fclaw_opt;
   fclaw_clawpatch_options_t *clawpatch_opt;
    fc3d_clawpack46_options_t   *claw46_opt;

    fclaw_global_t         *glob;
    fclaw_domain_t         *domain;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    fclaw_opt                  = fclaw_options_register(app, NULL, "fclaw_options.ini");
    clawpatch_opt  = fclaw_clawpatch_3d_options_register(app, "clawpatch", "fclaw_options.ini");
    claw46_opt       = fc3d_clawpack46_options_register(app, "claw3", "fclaw_options.ini");
    user_opt                 = latlong_options_register(app, "fclaw_options.ini");

    /* Read configuration file(s) */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt, user_opt);

        glob = fclaw_global_new();
        fclaw_global_store_domain(glob, domain);

        fclaw_options_store            (glob, fclaw_opt);
        fclaw_clawpatch_options_store (glob, clawpatch_opt);
        fc3d_clawpack46_options_store    (glob, claw46_opt);
        latlong_options_store            (glob, user_opt);

        run_program(glob);
        
        fclaw_global_destroy (glob);
    }

    fclaw_app_destroy (app);

    return 0;
}



