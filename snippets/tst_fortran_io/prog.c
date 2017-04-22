#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#define  MASTER     0

void print_hello_(int* mpirank);

int main (int argc, char *argv[])
{
    int  numtasks, mpirank, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
    MPI_Get_processor_name(hostname, &len);

    print_hello_(&mpirank);

    // printf ("Hello from task %d on %s!\n", taskid, hostname);
    if (mpirank == MASTER)
    {
        printf("MASTER: Number of MPI tasks is: %d\n",numtasks);
    }
    MPI_Finalize();
}

