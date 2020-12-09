#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <mpi.h>

int main(int argc, char *argv[])
{
	// global setup...
	MPI_Init(nullptr, nullptr);

	int result = Catch::Session().run(argc, argv);

	// abort if failure, some tests can hang otherwise
	if (result > 0) {
		MPI_Abort(MPI_COMM_WORLD, result);
	}

	// global clean-up...
	MPI_Finalize();

	return result;
}
