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

#include "overpressure_user.h"

#include "../all/euler3d_user.h"

static
void create_domain(fclaw_global_t *glob)
{
    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    int mi,mj,a,b;
    mi = fclaw_opt->mi;
    mj = fclaw_opt->mj;
    a = fclaw_opt->periodic_x;
    b = fclaw_opt->periodic_y;

    /* Used locally */
    double pi = M_PI;
    double rotate[2];

    rotate[0] = pi*fclaw_opt->theta/180.0;
    rotate[1] = pi*fclaw_opt->phi/180.0;

    /* Mapped, multi-block domain */
    fclaw_domain_t         *domain;
    fclaw_map_context_t    *cont = NULL, *brick = NULL;

    fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    fc3d_clawpack46_options_t *claw3_opt = fc3d_clawpack46_get_options(glob);

    /* Use [ax,bx]x[ay,by] */
    const user_options_t *user = overpressure_get_options(glob);
    switch (user->mapping) 
    {
    case 0:
        /* Use [ax,bx]x[ay,by]x[az,bz] to specify domain size */
        FCLAW_ASSERT(claw3_opt->mcapa == 0);
        FCLAW_ASSERT(fclaw_opt->manifold == 0);
        domain =
            fclaw_domain_new_2d_brick (glob->mpicomm, mi, mj, a, b,
                                     fclaw_opt->minlevel);

        /* Size is set by [ax,bx] x [ay, by], set in .ini file */            
        brick = fclaw_map_new_2d_brick(domain,mi,mj,a,b);        
        cont = fclaw_map_new_nomap_brick(brick);

        break;
    case 1:
        /* Square brick domain */
        FCLAW_ASSERT(claw3_opt->mcapa != 0);
        FCLAW_ASSERT(fclaw_opt->manifold != 0);
        FCLAW_ASSERT(clawpatch_opt->maux == 31);

        domain =
            fclaw_domain_new_2d_brick (glob->mpicomm, mi, mj, a, b,
                                     fclaw_opt->minlevel);

        brick = fclaw_map_new_2d_brick(domain,mi,mj,a,b);

        /* Square in [-1,1]x[-1,1], scaled/shifted to [0,1]x[0,1] */
        cont = fclaw2d_map_new_cart(brick,
                                    fclaw_opt->scale,
                                    fclaw_opt->shift);

        break;
    case 2:
        /* Latlong */
        /* Create brick domain with periodicity */        
        domain = fclaw_domain_new_2d_brick(glob->mpicomm, mi, mj, a, b,
                                          fclaw_opt->minlevel);

        /* Create brick mapping */        
        brick = fclaw_map_new_2d_brick(domain, mi, mj, a, b);
        
        cont = fclaw2d_map_new_latlong(brick,fclaw_opt->scale,
                                       rotate,
                                       user->latitude, 
                                       user->longitude,
                                       a,b);

        break;
    case 3:
        domain =
            fclaw_domain_new_2d_cubedsphere (glob->mpicomm,
                                            fclaw_opt->minlevel);
        cont = fclaw2d_map_new_cubedsphere (fclaw_opt->scale, rotate);

        break;
    case 4:
        domain =
            fclaw_domain_new_2d_twosphere (glob->mpicomm, fclaw_opt->minlevel);
        cont = fclaw2d_map_new_pillowsphere (fclaw_opt->scale, rotate);
        break;
    default:
        SC_ABORT_NOT_REACHED ();
    }

    if (user->mapping > 0) 
    {
        overpressure_map_extrude(cont,user->maxelev,
                                 user->mapping,
                                 user->min_z,
                                 user->max_z,
                                 user->mid_z,
                                 user->scale_bump);        
    }


    /* Store mapping in the glob */
    fclaw_map_store (glob, cont);            

    /* Store the domain in the glob */
    fclaw_global_store_domain(glob, domain);

    /* print out some info */    
    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
}

static
void run_program(fclaw_global_t* glob)
{
    const user_options_t *user_opt = overpressure_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (user_opt->claw_version == 4)
    {
        fc3d_clawpack46_solver_initialize(glob);
    }
    else if (user_opt->claw_version == 5)
    {
        printf("overpressure.cpp : Example not implemented for Claw version 5.\n");
        exit(0);
    }

    overpressure_link_solvers(glob);

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
    /* Initialize application */
    fclaw_app_t *app = fclaw_app_new (&argc, &argv, NULL);;

    /* Options */
    fclaw_options_t               *fclaw_opt;
    fclaw_clawpatch_options_t  *clawpatch_opt;
    fc3d_clawpack46_options_t     *claw46_opt;
    user_options_t                *user_opt;

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt =   fclaw_clawpatch_3d_options_register(app, "clawpatch",  "fclaw_options.ini");
    claw46_opt =        fc3d_clawpack46_options_register(app, "claw3", "fclaw_options.ini");
    user_opt =              overpressure_options_register(app,"fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    int first_arg;
    fclaw_exit_type_t vexit = 
        fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        /* Options have been checked and are valid */
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);
        fclaw_global_t *glob = fclaw_global_new_comm (mpicomm, size, rank);

        /* Store option packages in glob */
        fclaw_options_store           (glob, fclaw_opt);
        fclaw_clawpatch_options_store(glob, clawpatch_opt);
        fc3d_clawpack46_options_store   (glob, claw46_opt);
        overpressure_options_store       (glob, user_opt);

        /* Create domain and store domain in glob */
        create_domain(glob);

        /* Run the program */
        run_program(glob);
        
        fclaw_global_destroy(glob);
    }
    
    fclaw_app_destroy (app);

    return 0;
}



