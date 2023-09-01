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

#include "filament_user.h"

#include "../all/advection_user.h"

static
fclaw_domain_t* create_domain(sc_MPI_Comm mpicomm, 
                                fclaw_options_t* fclaw_opt, 
                                user_options_t* user,
                               fclaw_clawpatch_options_t* clawpatch_opt,
                                fc3d_clawpack46_options_t *claw3_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
    
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;    
    int a = 0; /* non-periodic */
    int b = 0;

    int mx = clawpatch_opt->d3->mx;
    int minlevel = fclaw_opt->minlevel;
    int check = mi*mx*pow_int(2,minlevel);

    switch (user->example) 
    {
    case 0:
        FCLAW_ASSERT(claw3_opt->mcapa == 0);
        FCLAW_ASSERT(fclaw_opt->manifold == 0);

        /* Size is set by [ax,bx] x [ay, by], set in .ini file */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        cont = fclaw2d_map_new_nomap_brick(brick);
        break;

    case 1:
        /* Square brick domain */
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
        /* Square brick domain */
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        /* Square in [-1,1]x[-1,1], shifted by (1,1,0) */
        cont = fclaw2d_map_new_cart(brick,
                                    fclaw_opt->scale,
                                    fclaw_opt->shift);

        break;
    case 2:
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
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
        FCLAW_ASSERT(mi == 2 && mj == 2);
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
        conn = p4est_connectivity_new_brick(mi,mj,a,b);
        brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
        cont = fclaw2d_map_new_bilinear (brick, 
                                         fclaw_opt->scale,
                                         fclaw_opt->shift, 
                                         user->center);
        break;

    default:
        SC_ABORT_NOT_REACHED ();
    }

    if (user->example > 0)
    {
        filament_map_extrude(cont,user->maxelev);
    }
    
    domain = fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont);
    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
    
    return domain;    
}

static
void run_program(fclaw_global_t* glob)
{
    user_options_t             *user;

    user = (user_options_t*) filament_get_options(glob);


    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw_domain_data_new(glob->domain);

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    if (user->claw_version == 4)
    {
      fc3d_clawpack46_solver_initialize(glob);
    }
    else if (user->claw_version == 5)
    {
        printf("filament.cpp : 3d example not implemented for Claw version 5.\n");
        exit(0);
    }

    filament_link_solvers(glob);

    fclaw_initialize(glob);
    fclaw_run(glob);
    fclaw_finalize(glob);
}



int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    user_options_t               *user_opt;
    fclaw_options_t              *fclaw_opt;
   fclaw_clawpatch_options_t *clawpatch_opt;
    fc3d_clawpack46_options_t    *claw46_opt;

    fclaw_global_t         *glob;
    fclaw_domain_t         *domain;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    fclaw_opt                    = fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt   = fclaw_clawpatch_options_register_3d(app, "clawpatch",  "fclaw_options.ini");
    claw46_opt         = fc3d_clawpack46_options_register(app, "claw3",      "fclaw_options.ini");
    user_opt =                  filament_options_register(app,               "fclaw_options.ini");  

    /* Read configuration file(s) */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, 
                               fclaw_opt, 
                               user_opt,
                               clawpatch_opt,
                               claw46_opt);
            
        glob = fclaw_global_new();
        fclaw_global_store_domain(glob, domain);

        fclaw_options_store            (glob, fclaw_opt);
        fclaw_clawpatch_options_store  (glob, clawpatch_opt);
        fc3d_clawpack46_options_store    (glob, claw46_opt);
        filament_options_store           (glob, user_opt);

        run_program(glob);

        fclaw_global_destroy(glob);        
    }

    fclaw_app_destroy (app);

    return 0;
}