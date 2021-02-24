#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "catch_reporter_automake.hpp"
#include <fclaw_mpi.h>

int main(int argc, char *argv[])
{
	// global setup...
	fclaw_mpi_init(nullptr, nullptr, MPI_COMM_WORLD, SC_LP_PRODUCTION);

	int result = Catch::Session().run(argc, argv);

	// abort if failure, some tests can hang otherwise
	if (result > 0) {
		return -1;
	}

	// global clean-up...
	fclaw_mpi_finalize();

	return result;
}
