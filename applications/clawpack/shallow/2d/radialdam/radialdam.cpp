/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#include "radialdam_user.h"

static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, 
                                fclaw_options_t* fclaw_opt,
                                user_options_t* user)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL;

    /* Local variables */
    double rotate[2];

    rotate[0] = fclaw_opt->phi;
    rotate[1] = fclaw_opt->theta;

    switch (user->example)
    {
    case 0:
        /* Use [ax,bx]x[ay,by] */
        conn = p4est_connectivity_new_unitsquare();
        cont = fclaw2d_map_new_nomap();
        break;
    case 1:
        /* Five patch square : maps to [-1,1]x[-1,1] */
        conn = p4est_connectivity_new_disk (0,0);
        cont = fclaw2d_map_new_fivepatch (fclaw_opt->scale,
                                          fclaw_opt->shift, user->alpha);
        break;
    case 2:
        /* Pillow disk (single block) */
        conn = p4est_connectivity_new_unitsquare();
        cont = fclaw2d_map_new_pillowdisk (fclaw_opt->scale,fclaw_opt->shift,rotate);
        break;
    case 3:
        /* Pillow disk mapping of the five patch square */
        conn = p4est_connectivity_new_disk(0, 0);
        cont = fclaw2d_map_new_pillowdisk5(fclaw_opt->scale,fclaw_opt->shift,
                                          rotate,user->alpha);
        break;
    default:
        SC_ABORT_NOT_REACHED ();
    }

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

    user_opt = radialdam_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize();
    }
    else if (user_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize();
    }

    radialdam_link_solvers(glob);
    
    fclaw_fpe_handling(glob);

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
    user_opt =                radialdam_options_register(app,"fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    radialdam_global_post_process(fclaw_opt, clawpatch_opt, user_opt);
    fclaw_app_print_options(app);

    /* Run the program */
    if (!retval & (vexit < 2))
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt,user_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw2d_options_store           (glob, fclaw_opt);
        fclaw2d_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_clawpack46_options_store   (glob, claw46_opt);
        fc2d_clawpack5_options_store    (glob, claw5_opt);
        radialdam_options_store         (glob, user_opt);

        run_program(glob);
        
        fclaw2d_global_destroy(glob);
    }
    
    fclaw_app_destroy (app);

    return 0;
}

