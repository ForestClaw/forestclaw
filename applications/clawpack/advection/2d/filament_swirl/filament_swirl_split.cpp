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

#include "filament/filament_user.h"
#include "swirl/swirl_user.h"


static
void run_programs(fclaw2d_global_t* globs[], int nglobs)
{
    for(int i = 0; i < nglobs; i++)
    {
        if(globs[i]->mpicomm != sc_MPI_COMM_NULL)
        {
            fclaw2d_set_global_context(globs[i]);

            fclaw2d_problem_setup(globs[i]);

            fclaw2d_clear_global_context(globs[i]);
        }
    }

    // Check if communicators are split
    int num_comms = 0; /* number of comms on this rank*/
    int glob_index = -1; /* if split, index of glob on this rank, otherwise value is garbage */
    for(int i = 0; i < nglobs; i++)
    {
        if(globs[i]->mpicomm != sc_MPI_COMM_NULL){
            num_comms += 1;
            glob_index = i;
        }
    }
    int one_comm = (num_comms == 1);
    int split_comms;
    sc_MPI_Allreduce(&one_comm, &split_comms, 1, \
                     sc_MPI_INT, sc_MPI_MIN, sc_MPI_COMM_WORLD);

    if (split_comms)
    {
        fclaw2d_set_global_context(globs[glob_index]);

        fclaw2d_run(globs[glob_index]);

        fclaw2d_clear_global_context(globs[glob_index]);
    }
    else
    {
        for(int i = 0; i < nglobs; i++)
        {
            fclaw2d_set_global_context(globs[i]);

            fclaw2d_run(globs[i]);

            fclaw2d_clear_global_context(globs[i]);
        }
    }
}

int
main (int argc, char **argv)
{
    /* Options */

    sc_MPI_Init(&argc, &argv);

    /* Split MPI_COMM_WORLD */
    sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

    int global_rank;
    sc_MPI_Comm_rank(mpicomm, &global_rank);

    int global_size;
    sc_MPI_Comm_size(mpicomm, &global_size);

    int color = global_rank < global_size/2 ? 0 : 1;
#if 0
    if(global_rank < global_size/2)
    {
        color = 0;
    }
    else
    {
        color = 1;
    }
#endif    

    sc_MPI_Comm subcomm;
    sc_MPI_Comm_split(mpicomm, color, global_rank, &subcomm);


    /* Initialize application on subcommunicator */
    fclaw_app_t *app = fclaw_app_new_on_comm(subcomm, &argc, &argv, NULL);


    /* Register packages */
    filament_options_t          *filament_user_opt;
    fclaw_options_t             *filament_fclaw_opt;
    fclaw2d_clawpatch_options_t *filament_clawpatch_opt;
    fc2d_clawpack46_options_t   *filament_claw46_opt;
    fc2d_clawpack5_options_t    *filament_claw5_opt;

    filament_fclaw_opt                    = 
        fclaw_options_register(app,  "filament",           "fclaw_options.ini");
    filament_clawpatch_opt    = 
        fclaw2d_clawpatch_options_register(app, "filament-clawpatch",  "fclaw_options.ini");
    filament_claw46_opt         = 
        fc2d_clawpack46_options_register(app, "filament-clawpack46", "fclaw_options.ini");
    filament_claw5_opt           = 
        fc2d_clawpack5_options_register(app, "filament-clawpack5",  "fclaw_options.ini");
    filament_user_opt =                  
        filament_options_register(app, "filament-user",       "fclaw_options.ini");  


    swirl_options_t             *swirl_user_opt;
    fclaw_options_t             *swirl_fclaw_opt;
    fclaw2d_clawpatch_options_t *swirl_clawpatch_opt;
    fc2d_clawpack46_options_t   *swirl_claw46_opt;
    fc2d_clawpack5_options_t    *swirl_claw5_opt;

    swirl_fclaw_opt =                   
        fclaw_options_register(app, "swirl",            "fclaw_options.ini");
    swirl_clawpatch_opt =   
        fclaw2d_clawpatch_options_register(app, "swirl-clawpatch",  "fclaw_options.ini");
    swirl_claw46_opt =        
        fc2d_clawpack46_options_register(app, "swirl-clawpack46", "fclaw_options.ini");
    swirl_claw5_opt =          
        fc2d_clawpack5_options_register(app, "swirl-clawpack5",  "fclaw_options.ini");
    swirl_user_opt =                    
        swirl_options_register(app, "swirl-user",       "fclaw_options.ini");  

    /* Read configuration file(s) */
    int first_arg;
    fclaw_exit_type_t vexit = 
        fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!vexit)
    {
        /* Options have been checked and are valid */

        /* MPI COMMs */
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);

        /* Globs */
        fclaw2d_global_t *filament_glob = fclaw2d_global_new_comm (subcomm, size, rank);

        fclaw2d_options_store            (filament_glob, filament_fclaw_opt);
        fclaw2d_clawpatch_options_store  (filament_glob, filament_clawpatch_opt);
        fc2d_clawpack46_options_store    (filament_glob, filament_claw46_opt);
        fc2d_clawpack5_options_store     (filament_glob, filament_claw5_opt);
        filament_options_store           (filament_glob, filament_user_opt);


        fclaw2d_global_t *swirl_glob = fclaw2d_global_new_comm (subcomm, size, rank);

        fclaw2d_options_store           (swirl_glob, swirl_fclaw_opt);
        fclaw2d_clawpatch_options_store (swirl_glob, swirl_clawpatch_opt);
        fc2d_clawpack46_options_store   (swirl_glob, swirl_claw46_opt);
        fc2d_clawpack5_options_store    (swirl_glob, swirl_claw5_opt);
        swirl_options_store             (swirl_glob, swirl_user_opt);

        if (color == 0)
        {
            /* Domain */
            // filament_domain = filament_create_domain(subcomm, filament_fclaw_opt, filament_user_opt,filament_clawpatch_opt);
            // fclaw2d_global_store_domain(filament_glob, filament_domain);

            filament_create_domain(filament_glob);

            /* initialize */
            filament_initialize(filament_glob);

            /* run */
            run_programs(&filament_glob, 1);

            /* finalize */
            filament_finalize(filament_glob);
        }
        else
        {
            /* Domain */
            //swirl_domain = swirl_create_domain(subcomm, swirl_fclaw_opt);
            //fclaw2d_global_store_domain(swirl_glob, swirl_domain);

            swirl_create_domain(swirl_glob);

            /* initialize */
            swirl_initialize(swirl_glob);

            /* run */
            run_programs(&swirl_glob, 1);

            /* finalize */
            swirl_finalize(swirl_glob);
        }

        /* destroy */
        fclaw2d_global_destroy(filament_glob);        
        fclaw2d_global_destroy(swirl_glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
