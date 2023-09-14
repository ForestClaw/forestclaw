/*
  Copyright (c) 2019-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "poisson_user.h"
    
#include <fclaw_include_all.h>

#include <fclaw_output.h>
#include <fclaw_diagnostics.h>

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.h>

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_options.h>


static
fclaw_domain_t* create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
 
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;

    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    /* Map unit square to disk using mapc2m_disk.f */
    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
    cont = fclaw2d_map_new_nomap_brick(brick);

    domain = fclaw_domain_wrap_2d(fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont));
    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
    return domain;
}

static
void run_program(fclaw_global_t* glob)
{
    // const poisson_options_t           *user_opt;

    // user_opt = poisson_get_options(glob);

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    /* Test thunderegg solver */
    fc2d_thunderegg_solver_initialize(glob);

    /* set up elliptic solver to use the thunderegg solver */
    poisson_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */

    /* Set up grid and RHS */
    fclaw_initialize(glob);

    /* Compute sum of RHS; reset error accumulators */
    int init_flag = 1;  
    fclaw_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    /* Output rhs */
    int Frame = 0;
    fclaw_output_frame(glob,Frame);
 
    /* Solve the elliptic problem */
    fclaw_elliptic_solve(glob);

    /* Compute error, compute conservation */
    fclaw_diagnostics_gather(glob, init_flag);                

    /* Output solution */
    Frame = 1;
    fclaw_output_frame(glob,Frame);

    /* ---------------------------------------------------------------
       Finalize
       --------------------------------------------------------------- */
    fclaw_finalize(glob);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    fclaw_options_t             *fclaw_opt;

    fclaw_clawpatch_options_t *clawpatch_opt;
    fc2d_thunderegg_options_t    *mg_opt;
    poisson_options_t              *user_opt;

    fclaw_global_t            *glob;
    fclaw_domain_t            *domain;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt =   fclaw_clawpatch_2d_options_register(app, "clawpatch",  "fclaw_options.ini");
    mg_opt =            fc2d_thunderegg_options_register(app, "thunderegg", "fclaw_options.ini");
    user_opt =                  poisson_options_register(app,               "fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw_global_new();
        fclaw_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw_options_store           (glob, fclaw_opt);
        fclaw_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_thunderegg_options_store    (glob, mg_opt);
        poisson_options_store            (glob, user_opt);

        run_program(glob);

        fclaw_global_destroy(glob);        
    }
    
    fclaw_app_destroy (app);

    return 0;
}
