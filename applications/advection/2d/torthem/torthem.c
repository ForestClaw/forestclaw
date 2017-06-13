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

#include "torus_common.h"

#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack46.h>
#include <fclaw2d_clawpatch.h>


typedef struct torthem
{
    fclaw2d_global_t *global;

    fc2d_clawpack46_options_t claw_opt;
    fclaw2d_clawpatch_options_t clawpatch_opt;
    fclaw_options_t fclaw_opt;
    user_options_t user;
}
torthem_t;

static void
torthem_init (torthem_t * torthem)
{
    fclaw_options_t *fclaw_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;    
    fc2d_clawpack46_options_t *claw_opt;
    user_options_t *user;
    fclaw_exit_type_t et;

    memset (torthem, 0, sizeof (*torthem));
    torthem->global = fclaw2d_global_new (NULL);

    claw_opt      = &torthem->claw_opt;
    user          = &torthem->user;
    fclaw_opt     = &torthem->fclaw_opt;
    clawpatch_opt = &torthem->clawpatch_opt;

    /**************** FCLAW OPTIONS *************/

    /* todo: set all other values */
    fclaw_opt->dim = 2;
    fclaw_opt->periodic_x = 1;
    fclaw_opt->periodic_y = 0;
    fclaw_opt->mi = 18;
    fclaw_opt->mj = 5;
    fclaw_opt->manifold = 1;

    fclaw_opt->outstyle = 3;
    fclaw_opt->nout = 10;
    fclaw_opt->nstep = 1;

    fclaw_opt->trapfpe = 0;
    fclaw_opt->mpi_debug = 0;

    fclaw_opt->verbosity = 0;

    fclaw_opt->scale_string = "1. 1. 1.";
    fclaw_opt->shift_string = "0. 0. 0.";

    et = fclaw_options_postprocess (fclaw_opt);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Option postprocess error");
    et = fclaw_options_check (fclaw_opt);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Option check error");

    /**************** CLAWPACK 4.6 OPTIONS *************/
    clawpatch_opt->mx = 8;
    clawpatch_opt->my = 8;
    clawpatch_opt->mbc = 2;
    clawpatch_opt->meqn = 1;

    /**************** CLAWPACK 4.6 OPTIONS *************/

    claw_opt->mthbc_string = "1 1 1 1";

    claw_opt->mwaves = 3;
    claw_opt->mthlim_string = "1 1 4";
    claw_opt->order_string = "2 2";
    /* Plus more options */

    et = fc2d_clawpack46_postprocess (claw_opt);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Clawpack46 postprocess error");

    et = fc2d_clawpack46_check (claw_opt);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Clawpack46 check error");

    for(int i = 0; i < claw_opt->mwaves; i++)
    {
        fclaw_global_essentialf("mthlim[%d] = %d\n",i,claw_opt->mthlim[i]);
    }

    fc2d_clawpack46_vtable_initialize();

    /**************** TORUS OPTIONS *************/

    user->latitude_string = "-50. 50.";
    user->longitude_string = "0. 360.";

    et = torus_options_postprocess (user);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Torus postprocess error");
    et = torus_options_check (user);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Torus check error");
}

static void
torthem_run (torthem_t * torthem)
{
    fclaw_global_essentialf("Running torthem ...\n");
    fclaw_global_essentialf("Done.\n");
}

static void
torthem_destroy (torthem_t * torthem)
{
    torus_options_destroy (&torthem->user);
    fc2d_clawpack46_options_destroy (&torthem->claw_opt);
    fclaw_options_destroy (torthem->fclaw_opt);
    fclaw2d_global_destroy (torthem->global);
}

int
main (int argc, char **argv)
{
    int mpiret;
    sc_MPI_Comm mpicomm;

    torthem_t storthem, *torthem = &storthem;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init (NULL, SC_LP_ESSENTIAL);
    fclaw_init (NULL, SC_LP_INFO);

    torthem_init (torthem);
    torthem_run (torthem);
    torthem_destroy (torthem);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}
