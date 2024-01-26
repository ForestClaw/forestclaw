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

/* This example demonstrates the use of fclaw2d_overlap_exchange to exchange
 * interpolation data between two meshes.
 * Both meshes have a clearly assigned role:
 *  consumer (e.g. Gemini) - queries data for points - represented by swirl
 *  producer (e.g. MAGIC) - provides data - represented by filament. */

#include "filament/filament_user.h"
#include "swirl/swirl_user.h"

#include "user.h"

#include "overlap.h"

#include <fclaw_global.h>

static
void setup_overlap(fclaw_global_t *swirl_glob,\
                   fclaw_global_t *filament_glob,
                   overlap_consumer_t *c)
{
    /* -------------- setup overlap infomation -------------------*/
    /* compute process-local query points on the consumer side */

    fclaw2d_clawpatch_options_t *swirl_clawpatch_opt = 
        fclaw2d_clawpatch_get_options(swirl_glob);

    //overlap_consumer_t consumer, *c = &consumer;
    c->glob = swirl_glob;
    c->domain = swirl_glob->domain;
    int mx = swirl_clawpatch_opt->mx;
    int my = swirl_clawpatch_opt->my;
    c->num_cells_in_patch = mx*my;
    create_query_points (c);

    /* initialize the filament geometry information that is needed for
     * mapping between the swirl and the filament domain */
    fclaw_options_t *filament_fclaw_opt = fclaw_get_options(filament_glob);
    overlap_geometry_t filament_geometry, *geo = &filament_geometry;

    geo->fclaw_opt = filament_fclaw_opt;
    geo->blocks = filament_glob->domain->blocks;

    /* obtain interpolation data of the points from the producer side */
    fclaw_overlap_exchange (filament_glob->domain, c->query_points,
                              overlap_interpolate, geo);

    /* output the interpolation data for all query points */
    output_query_points (c);
}

static
void run_program(fclaw_global_t *swirl_glob,\
                fclaw_global_t *filament_glob)
{
    /* run */
    fclaw_global_t *globs[2];
    globs[0] = filament_glob;
    globs[1] = swirl_glob;
    user_run (globs, 2);

    /* Finalize solvers */
    filament_finalize(filament_glob);
    swirl_finalize(swirl_glob);    
}


int
main (int argc, char **argv)
{
    /* Initialize application */
    fclaw_app_t *app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    //Global options like verbosity, etc
    fclaw_app_options_register_core(app, "filament_options.ini"); 

    /* Filament options */
    filament_options_t          *filament_user_opt;
    fclaw_options_t             *filament_fclaw_opt;
    fclaw2d_clawpatch_options_t *filament_clawpatch_opt;
    fc2d_clawpack46_options_t   *filament_claw46_opt;
    fc2d_clawpack5_options_t    *filament_claw5_opt;

    filament_fclaw_opt = 
                    fclaw_options_register(app, "filament",            "filament_options.ini");
    filament_clawpatch_opt    = 
        fclaw2d_clawpatch_options_register(app, "filament-clawpatch",  "filament_options.ini");
    filament_claw46_opt = 
          fc2d_clawpack46_options_register(app, "filament-clawpack46", "filament_options.ini");
    filament_claw5_opt           = 
           fc2d_clawpack5_options_register(app, "filament-clawpack5",  "filament_options.ini");
    filament_user_opt =                  
                 filament_options_register(app, "filament-user",       "filament_options.ini");  

    /* Swirl options */
    swirl_options_t             *swirl_user_opt;
    fclaw_options_t             *swirl_fclaw_opt;
    fclaw2d_clawpatch_options_t *swirl_clawpatch_opt;
    fc2d_clawpack46_options_t   *swirl_claw46_opt;
    fc2d_clawpack5_options_t    *swirl_claw5_opt;

    swirl_fclaw_opt =                   
                    fclaw_options_register(app, "swirl",           "swirl_options.ini");
    swirl_clawpatch_opt =   
        fclaw2d_clawpatch_options_register(app, "swirl-clawpatch", "swirl_options.ini");
    swirl_claw46_opt =        
          fc2d_clawpack46_options_register(app, "swirl-clawpack46","swirl_options.ini");
    swirl_claw5_opt =          
           fc2d_clawpack5_options_register(app, "swirl-clawpack5",  "swirl_options.ini");
    swirl_user_opt =                    
                    swirl_options_register(app, "swirl-user",       "swirl_options.ini");  

    /* Read configuration file(s) */
    int first_arg;
    fclaw_exit_type_t vexit = 
        fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        /* Options have been checked and are valid */
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);

        /* Filament setup */
        fclaw_global_t *filament_glob = fclaw_global_new_comm (mpicomm, size, rank);

        fclaw_options_store            (filament_glob, filament_fclaw_opt);
        fclaw2d_clawpatch_options_store  (filament_glob, filament_clawpatch_opt);
        fc2d_clawpack46_options_store    (filament_glob, filament_claw46_opt);
        fc2d_clawpack5_options_store     (filament_glob, filament_claw5_opt);
        filament_options_store           (filament_glob, filament_user_opt);

        filament_create_domain(filament_glob);

        /* Swirl setup */
        fclaw_global_t *swirl_glob = fclaw_global_new_comm (mpicomm, size, rank);
            
        fclaw_options_store           (swirl_glob, swirl_fclaw_opt);
        fclaw2d_clawpatch_options_store (swirl_glob, swirl_clawpatch_opt);
        fc2d_clawpack46_options_store   (swirl_glob, swirl_claw46_opt);
        fc2d_clawpack5_options_store    (swirl_glob, swirl_claw5_opt);
        swirl_options_store             (swirl_glob, swirl_user_opt);

        swirl_create_domain(swirl_glob);

        /* initialize both solvers before doing overlap */
        filament_initialize(filament_glob);
        swirl_initialize(swirl_glob);    

        /* Set up overlap points */
        overlap_consumer_t consumer, *c = &consumer;
        setup_overlap(swirl_glob, filament_glob,c);

        /* Run the program */
        run_program(swirl_glob, filament_glob);

        /* destroy overlap points */
        sc_array_destroy (c->query_points);

        /* Destroy apps */
        fclaw_global_destroy (filament_glob);
        fclaw_global_destroy (swirl_glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
