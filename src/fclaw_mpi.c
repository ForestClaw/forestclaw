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

#include <fclaw_mpi.h>
#include <fclaw2d_global.h>

/* Functions with C prototypes to use forestclaw from C code */

#ifndef P4EST_ENABLE_MPI
static int initialized = 0;
#endif

void
fclaw_mpi_init (int * argc, char *** argv, sc_MPI_Comm mpicomm, int lp)
{
#ifdef P4EST_ENABLE_MPI
    int mpiret;

    //mpiret = sc_MPI_Init (argc, argv);
    //SC_CHECK_MPI (mpiret);

    int provided;
    mpiret = sc_MPI_Init_thread (argc, argv, sc_MPI_THREAD_FUNNELED, &provided);
    if (provided != sc_MPI_THREAD_FUNNELED) {
        printf("Recieved mpi_init_thread level %d\n", provided);
    }
    SC_CHECK_MPI (mpiret);
#else
    initialized = 1;
#endif

    sc_init (mpicomm, 0, 0, NULL, lp);
    p4est_init (NULL, lp);
}

void
fclaw_mpi_initialized (int * flag)
{
#ifdef P4EST_ENABLE_MPI
    MPI_Initialized(flag);
#else
    *flag = initialized;
#endif
}

void
fclaw_mpi_finalize (void)
{
    int mpiret;

    sc_finalize ();

    /* With P4EST_ENABLE_MPI, MPI_Init and MPI_Finalize must match up.
     * Without it does not matter. */
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
}

void
fclaw_mpi_barrier (sc_MPI_Comm mpicomm)
{
    int mpiret;

    mpiret = sc_MPI_Barrier (mpicomm);
    SC_CHECK_MPI (mpiret);
}

void
fclaw_mpi_serialization_enter (struct fclaw2d_global *glob)
{
    fclaw2d_domain_serialization_enter(glob->domain);
}

void
fclaw_mpi_serialization_leave (struct fclaw2d_global *glob)
{
    fclaw2d_domain_serialization_leave(glob->domain);
}



void
fclaw_mpi_debug (void)
{
#ifdef FCLAW_ENABLE_MPI
  int i;

  /* Find out process rank */
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out number of processes */
  int num;
  MPI_Comm_size(MPI_COMM_WORLD, &num);

  /* Don't do anything until we are done with this part */
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0)
  {
      printf("\n");
      printf("Getting setup for parallel debugging\n");
      printf("------------------------------------\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(i = 0; i < num; i++)
  {
      if (my_rank == i)
      {
          printf("Proc %d with process %d is waiting to be attached\n",my_rank,getpid());
          fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
  }

  int ii = 0;
  while (ii == 0)  /* (gdb) set ii=1 */
  {
    /* Code will stop here;  set ii=1 to continue in gdb */
  }

#else
    fclaw_global_infof
      ("Option --mpi_debug is set, but code is not compiled for MPI");
#endif
}
