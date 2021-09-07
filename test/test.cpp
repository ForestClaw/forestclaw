#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "catch_reporter_automake.hpp"
#include <fclaw_mpi.h>

int main(int argc, char *argv[])
{
	// global setup...
	fclaw_mpi_init(nullptr, nullptr, sc_MPI_COMM_WORLD, SC_LP_PRODUCTION);

	int result = Catch::Session().run(argc, argv);

	// abort if failure, some tests can hang otherwise
	if (result > 0) {
		return -1;
	}

	// global clean-up...
	//fclaw_mpi_finalize();

	return result;
}

#include <fclaw2d_domain.h>
#include <fclaw2d_options.h>
#include <fclaw2d_p4est.h>
#include <fclaw2d_convenience.h>

fclaw2d_domain_t* create_test_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
    
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;    
    int a = 0; /* non-periodic */
    int b = 0;

    /* Square brick domain */
    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick(conn,mi,mj);
    cont = fclaw2d_map_new_nomap_brick(brick);
    
    domain = fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
    
    return domain;    
}
