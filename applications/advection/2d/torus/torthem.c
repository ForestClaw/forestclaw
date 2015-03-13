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

static void
global_init (fclaw2d_global_t * glob)
{
    fclaw_options_t *gparms = &glob->gparms;
    fclaw_exit_type_t et;

    /* todo: set all other values */
    gparms->dim = 2;
    gparms->mx = 8;
    gparms->my = 8;
    gparms->mbc = 2;

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
}

static void
global_reset (fclaw2d_global_t * glob)
{
    fclaw_options_reset (&glob->gparms);
}

int
main (int argc, char **argv)
{
    int mpiret;
    sc_MPI_Comm mpicomm;

    fclaw2d_global_t global, *glob = &global;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init (NULL, SC_LP_ESSENTIAL);
    fclaw_init (NULL, SC_LP_INFO);

    global_init (glob);

    global_reset (glob);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}
