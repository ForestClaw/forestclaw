#include <atomic>
#include <fclaw_mpi.h>

/* Get whatever definitions exist already */
#ifdef FCLAW_HAVE_FENV_H
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>
#endif

/* Use as an alternate to GNU feenableexcept */
#ifndef FCLAW_HAVE_FEENABLEEXCEPT
#include <fp_exception_glibc_extension.h>
#endif

#define DOCTEST_CONFIG_IMPLEMENT
#include <test.hpp>
#include <exception>
#include <csetjmp>

static bool output_vtk=false;
static bool test_indirect_flag = false;

bool test_output_vtk()
{
    return output_vtk;
}

bool test_indirect()
{
    return test_indirect_flag;
}

static bool has_aborted=true;
static bool expect_abort=false;

std::jmp_buf jump_buffer;

void throw_exception()
{
    has_aborted = true;
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
    //add vtk option to output vtk files
    for (int i = 0; i < argc; i++) {
        output_vtk = strcmp(argv[i], "--vtk") == 0;
        if (output_vtk)
        {
            break;
        }
    }
    for (int i = 0; i < argc; i++) {
        test_indirect_flag = strcmp(argv[i], "--indirect") == 0;
        if (test_indirect_flag)
        {
            break;
        }
    }
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
        feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

        sc_set_abort_handler(throw_exception);
	    result = context.run();

	    // global clean-up...
	    sc_MPI_Finalize();

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