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

#include <fclaw2d_global.h>

#include "torus_common.h"

typedef struct torthem
{
    fclaw2d_global_t *global;
    user_options_t user;
}
torthem_t;

static void
torthem_init (torthem_t * torthem)
{
    fclaw_options_t *gparms;
    user_options_t *user = &torthem->user;
    fclaw_exit_type_t et;

    memset (torthem, 0, sizeof (*torthem));
    torthem->global = fclaw2d_global_new (NULL);
    gparms = torthem->global->gparms;

    /* todo: set all other values */
    gparms->dim = 2;
    gparms->mx = 8;
    gparms->my = 8;
    gparms->mbc = 2;

    gparms->periodic_x = 1;
    gparms->periodic_y = 0;
    gparms->mi = 18;
    gparms->mj = 5;
    gparms->manifold = 1;

    gparms->outstyle = 3;
    gparms->nout = 10;
    gparms->nstep = 1;

    gparms->trapfpe = 0;
    gparms->mpi_debug = 0;

    gparms->verbosity = 0;

    gparms->mthbc_string = "1 1 1 1";

    gparms->scale_string = "1. 1. 1.";
    gparms->shift_string = "0. 0. 0.";

    et = fclaw_options_postprocess (gparms);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Option postprocess error");
    et = fclaw_options_check (gparms);
    SC_CHECK_ABORT (et == FCLAW_NOEXIT, "Option check error");

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
}

static void
torthem_reset (torthem_t * torthem)
{
    torus_options_reset (&torthem->user);
    fclaw_options_reset (torthem->global->gparms);
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
    torthem_reset (torthem);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}
