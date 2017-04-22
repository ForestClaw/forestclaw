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
void print_hello_(int* mpirank, output_t* verb);

int main (int argc, char *argv[])
{
    int  numtasks, mpirank, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    output_t verb;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
    MPI_Get_processor_name(hostname, &len);

    verb = DEBUG;
    print_hello_1_(&mpirank,&verb);    // Doesn't really work

    // printf ("Hello from task %d on %s!\n", mpirank, hostname);
    if (mpirank == 0)
    {
        printf("Number of MPI tasks is: %d\n",numtasks);
    }
    MPI_Finalize();
}

