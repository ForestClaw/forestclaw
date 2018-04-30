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

#include "replicated_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack5_options.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

#include <math.h>

static int s_user_options_package_id = -1;

static void *
replicated_register (user_options_t *user_opt, sc_options_t * opt)
{
    sc_options_add_int (opt, 0, "example", &user_opt->example, 0,
                        "[user] 0 - multi-block; 1 - single block [0]");

    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 4,
                        "[user] Clawpack version (4 or 5) [5]");

    sc_options_add_int (opt, 0, "replicate-factor", &user_opt->replicate_factor, 1,
                        "[user] Replication factor [-1]");

    sc_options_add_int (opt, 0, "minlevel-base", &user_opt->minlevel_base, 4,
                        "[user] Block 1x1 minimum refinement [4]");

    sc_options_add_int (opt, 0, "maxlevel-base", &user_opt->maxlevel_base, 7,
                        "[user] Block 1x1 maximum refinement [7]");


    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t 
replicated_postprocess(user_options_t *user, fclaw_options_t *fclaw_opt)
{

    /* Check to make sure replicate factor is power of 2. 
       We do it here, since options_check only happens after post-process. */
    int pos = 0;    
    while (user->replicate_factor != (1 << pos))
    {
        pos++;
        if (pos > 32)
        {
            fclaw_global_essentialf("Replicated_postprocess : Replicate " \
                                    "factor should be a low power of 2 (< 2**32)\n");
            return FCLAW_EXIT_QUIET;
        }
    }

    fclaw_opt->ax = 0;
    fclaw_opt->ay = 0;
    fclaw_opt->bx = user->replicate_factor;
    fclaw_opt->by = user->replicate_factor;

    fclaw_opt->minlevel = user->minlevel_base;
    fclaw_opt->maxlevel = user->maxlevel_base;



    if (user->example == 0)
    {
        /* 1x1 block;  use replicate factor to dimension domain */
        fclaw_opt->mi = 1;
        fclaw_opt->mj = 1;

        /* We need to increase the refinement level when we increase domain size. 
           Code below computes pos=log2(replicate_factor) */
        fclaw_opt->minlevel += pos;
        fclaw_opt->maxlevel += pos;
    }
    else
    {
        /* NxN block */
        fclaw_opt->mi = user->replicate_factor;
        fclaw_opt->mj = user->replicate_factor;
    }

    fclaw_opt->periodic_x = 1;
    fclaw_opt->periodic_y = 1;
    fclaw_opt->smooth_level = fclaw_opt->maxlevel-1;
    return FCLAW_NOEXIT;
}

static void
replicated_destroy (user_options_t *user)
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

    return replicated_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;
    fclaw_options_t *fclaw_opt = 
                 (fclaw_options_t*) fclaw_app_get_attribute(app,"Options",NULL);

    return replicated_postprocess(user, fclaw_opt);
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

    replicated_destroy (user);

    FCLAW_FREE (user);
}


static const
fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    NULL,
    options_destroy
};

/* ------------- User options access functions --------------------- */

static
user_options_t* replicated_options_register (fclaw_app_t * app,
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

static 
void replicated_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

const user_options_t* replicated_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}
/* ------------------------- ... and here ---------------------------- */

static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, 
                                fclaw_options_t* fclaw_opt,
                                user_options_t* user_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

    int mi,mj,a,b;
    mi = fclaw_opt->mi;
    mj = fclaw_opt->mj;
    a = fclaw_opt->periodic_x;
    b = fclaw_opt->periodic_y;


    /* Duplicate initial conditions in each block */
    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick(conn,mi,mj);  /* this writes out brick data */
    cont = fclaw2d_map_new_nomap_brick(brick);

    domain = fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
    return domain;
}

static
void run_program(fclaw2d_global_t* glob)
{
    const user_options_t           *user_opt;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    user_opt = replicated_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtable_initialize();
    fclaw2d_diagnostics_vtable_initialize();

    /* Initialize virtual tables for solvers */
    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize();
    }
    else if (user_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize();
    }

    replicated_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
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
    sc_options_t                *options;
    user_options_t              *user_opt;
    fclaw_options_t             *fclaw_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;
    fc2d_clawpack46_options_t   *claw46_opt;
    fc2d_clawpack5_options_t    *claw5_opt;

    fclaw2d_global_t            *glob;
    fclaw2d_domain_t            *domain;
    sc_MPI_Comm mpicomm;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,"fclaw_options.ini");
    clawpatch_opt =   fclaw2d_clawpatch_options_register(app,"fclaw_options.ini");
    claw46_opt =        fc2d_clawpack46_options_register(app,"fclaw_options.ini");
    claw5_opt =          fc2d_clawpack5_options_register(app,"fclaw_options.ini");
    user_opt =               replicated_options_register(app,"fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!retval & !vexit)
    {

        /* Options have been checked and are valid */
        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt, user_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw2d_options_store           (glob, fclaw_opt);
        fclaw2d_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_clawpack46_options_store   (glob, claw46_opt);
        fc2d_clawpack5_options_store    (glob, claw5_opt);
        replicated_options_store        (glob, user_opt);

        /* Run the program */
        run_program(glob);
        
        fclaw2d_global_destroy(glob);
    }
    
    fclaw_app_destroy (app);

    return 0;
}


