#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

typedef enum
{
    SILENT=0,
    ESSENTIAL,
    PRODUCTION,
    INFO,
    DEBUG,
} output_t;


void print_hello_1_(int* mpirank, output_t* verb);
void print_hello_(int* mpirank);

int main (int argc, char *argv[])
{
    int  numprocs, mpirank, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    output_t verb;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Get_processor_name(hostname, &len);

    verb = DEBUG;
    print_hello_(&mpirank);

    if (mpirank == 0)
    {
        printf("Number of MPI tasks is: %d\n",numprocs);
    }
    MPI_Finalize();
}

