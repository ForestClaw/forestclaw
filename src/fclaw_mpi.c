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

#include <p4est_base.h>
#include "fclaw_mpi.h"

/* Functions with C prototypes to use forestclaw from C code */

void
fclaw_mpi_init (int * argc, char *** argv, sc_MPI_Comm mpicomm, int lp)
{
#ifdef P4EST_MPI
    int mpiret;

    //mpiret = sc_MPI_Init (argc, argv);
    //SC_CHECK_MPI (mpiret);

    int provided;
    mpiret = sc_MPI_Init_thread (argc, argv, sc_MPI_THREAD_FUNNELED, &provided);
    if (provided != sc_MPI_THREAD_FUNNELED) printf("Recieved mpi_init_thread level %d\n", provided);
    SC_CHECK_MPI (mpiret);

    sc_init (mpicomm, 0, 0, NULL, lp);
    p4est_init (NULL, lp);
#endif
}

void
fclaw_mpi_finalize (void)
{
    int mpiret;

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
}
