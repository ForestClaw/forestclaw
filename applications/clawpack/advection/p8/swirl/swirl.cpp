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

#include "swirl_user.h"

#include <fclaw3d_defs.h>

static
fclaw_domain_t* create_domain(sc_MPI_Comm mpicomm,
                                fclaw_options_t* fclaw_opt,
                                user_options_t *user,
                                void *clawpatch_opt,
                                void *claw3_opt)
{
    /* Mapped, multi-block domain */
    p8est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain = NULL;

#ifdef P8HACK
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
#endif

    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;
    int mk = 1;
    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;
    int c = 0;

    int minlevel = fclaw_opt->minlevel;

#ifdef P8HACK
    int mx = clawpatch_opt->mx;
    int check = mi*mx*pow_int(2,minlevel);
#endif

    switch (user->example)
    {
    case 0:
#ifdef P8HACK
        FCLAW_ASSERT(claw3_opt->mcapa == 0);
        FCLAW_ASSERT(fclaw_opt->manifold == 0);
#endif
        /* Size is set by [ax,bx] x [ay, by], set in .ini file */
        conn = p8est_connectivity_new_unitcube();
#ifdef P8HACK
        cont = fclaw2d_map_new_nomap();
#endif
        break;

    case 1:
        /* Square brick domain */
#ifdef P8HACK
        FCLAW_ASSERT(claw3_opt->mcapa != 0);
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
        FCLAW_ASSERT(clawpatch_opt->maux == 4);
#endif
        conn = p8est_connectivity_new_brick(mi,mj,mk,a,b,c);
#ifdef P8HACK
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        /* Square in [-1,1]x[-1,1], scaled/shifted to [0,1]x[0,1] */
        cont = fclaw2d_map_new_cart(brick,
                                    fclaw_opt->scale,
                                    fclaw_opt->shift);
#endif
        break;
#ifdef P8HACK
    case 2:
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
        FCLAW_ASSERT(claw3_opt->mcapa != 0);
        FCLAW_ASSERT(clawpatch_opt->maux == 4);
        if (check < 32)
        {
            printf("mi*pow_int(mx,minlevel) = %d\n",check);
            fclaw_global_essentialf("The five patch mapping requires mi*mx*2^minlevel > 32\n");
            exit(0);

        }
        /* Five patch square domain */
        conn = p4est_connectivity_new_disk (0, 0);
        cont = fclaw2d_map_new_fivepatch (fclaw_opt->scale,
                                          fclaw_opt->shift,
                                          user->alpha);

        break;

    case 3:
        /* bilinear square domain : maps to [-1,1]x[-1,1] */
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
        FCLAW_ASSERT(claw3_opt->mcapa != 0);
        FCLAW_ASSERT(clawpatch_opt->maux == 4);
        FCLAW_ASSERT(mi == 2 && mj == 2);
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        cont = fclaw2d_map_new_bilinear (brick,
                                         fclaw_opt->scale,
                                         fclaw_opt->shift,
                                         user->center);
        break;
#endif /* P8HACK */
    default:
        SC_ABORT_NOT_REACHED ();
    }

#ifdef P8HACK
    if (user->example > 0)
    {
        swirl_map_extrude(cont,user->maxelev);
    }
#endif /* P8HACK */

    domain = fclaw_domain_wrap_3d(fclaw3d_domain_new_conn (mpicomm, fclaw_opt->minlevel, conn));
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
    fclaw_domain_data_new(glob->domain);

    user_opt = swirl_get_options(glob);

#ifdef P8HACK
    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    if (user_opt->claw_version == 4)
    {
        fc3d_clawpack46_solver_initialize(glob);
    }
    else if (user_opt->claw_version == 5)
    {
        printf("swirl.cpp : Example not implemented for Claw version 5.\n");
        exit(0);
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw_initialize(glob);
    fclaw_run(glob);
    fclaw_finalize(glob);

#endif /* P8HACK */
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;

    /* Options */
    user_options_t               *user_opt;
#if 0
   fclaw_clawpatch_options_t *clawpatch_opt = NULL;
    fc3d_clawpack46_options_t    *claw46_opt = NULL;
#endif
    void                         *clawpatch_opt = NULL;
    void                         *claw46_opt = NULL;
    fclaw_options_t              *fclaw_opt;
    sc_options_t                 *options;

    int first_arg;
    fclaw_exit_type_t vexit;

    fclaw_global_t            *glob;
    fclaw_domain_t            *domain;

    int retval;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,  NULL,       "fclaw_options.ini");

#ifdef P8HACK
    clawpatch_opt =  fclaw_clawpatch_options_register_3d(app, "clawpatch", "fclaw_options.ini");
    claw46_opt =        fc3d_clawpack46_options_register(app, "claw3",     "fclaw_options.ini");
#endif /* P8HACK */

    user_opt =                    swirl_options_register(app,              "fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!retval & !vexit)
    {
        /* Options have been checked and are valid */
        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt, user_opt, clawpatch_opt,
                               claw46_opt);

        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw_global_new();
        fclaw_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw_options_store           (glob, fclaw_opt);

#ifdef P8HACK
        fclaw_clawpatch_options_store(glob, clawpatch_opt);
        fc3d_clawpack46_options_store   (glob, claw46_opt);
#endif /* P8HACK */
        swirl_options_store             (glob, user_opt);

        run_program(glob);

        fclaw_global_destroy(glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
