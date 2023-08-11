#include <atomic>
#include <fclaw_mpi.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <test.hpp>
#include <exception>
#include <csetjmp>

static bool expect_abort=false;

std::jmp_buf jump_buffer;

void throw_exception()
{
    if(expect_abort)
    {
        expect_abort=false;
        std::longjmp(jump_buffer, 1);
    }
}

void fclaw_test_expect_abort()
{
    expect_abort=true;
}

void fclaw_test_clear_expect_abort()
{
    expect_abort=false;
}

std::jmp_buf& fclaw_test_get_jump_buffer()
{
    return jump_buffer;
}

int main(int argc, char *argv[])
{
    bool listing = false;
    for (int i = 0; i < argc; i++) {
        listing = strcmp(argv[i], "--list-test-cases") == 0;
        if (listing)
            break;
        listing = strcmp(argv[i], "--list-test-suites") == 0;
        if (listing)
            break;
        listing = strcmp(argv[i], "--list-reporters") == 0;
        if (listing)
            break;
        listing = strcmp(argv[i], "-ltc") == 0;
        if (listing)
            break;
        listing = strcmp(argv[i], "-lts") == 0;
        if (listing)
            break;
        listing = strcmp(argv[i], "-ltr") == 0;
        if (listing)
            break;
    }

    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int result = 0;

    if(listing) {
        int rank = 0;
        if(rank == 0){
	       result = context.run();
        }
    } else {
	    // global setup...
	    fclaw_mpi_init(nullptr, nullptr, sc_MPI_COMM_WORLD, SC_LP_PRODUCTION);

        sc_set_abort_handler(throw_exception);
	    result = context.run();

	    // global clean-up...
	    //fclaw_mpi_finalize();

	    return result;
    } 

    // abort if failure, some tests can hang otherwise
    if (!listing && (result > 0 || context.shouldExit())) {
#ifdef P4EST_ENABLE_MPI
        MPI_Abort(MPI_COMM_WORLD, result);
#else
        return -1;
#endif
    }

    return 0;
}

#include <fclaw_domain.h>
#include <fclaw2d_options.h>
#include <fclaw2d_p4est.h>
#include <fclaw2d_convenience.h>

fclaw_domain_t* create_test_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
    
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;    
    int a = 0; /* non-periodic */
    int b = 0;

    /* Square brick domain */
    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick_conn (conn,mi,mj);
    cont = fclaw2d_map_new_nomap_brick(brick);
    
    domain = fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
    
    return domain;    
}
