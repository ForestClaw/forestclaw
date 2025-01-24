#include <atomic>
#include <fclaw_mpi.h>
#include <fclaw_forestclaw.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <test.hpp>
#include <exception>
#include <csetjmp>

static bool output_vtk=false;

bool test_output_vtk()
{
    return output_vtk;
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
    else{
        FAIL("An abort was called");
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

static bool dirty_memory=false;
void fclaw_test_set_dirty_memory()
{
    dirty_memory = true;
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
            std::cout << "outputting vtk files" << std::endl;
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
        fclaw_app_t * app = fclaw_app_new(&argc, &argv, nullptr);

        sc_set_abort_handler(throw_exception);
	    result = context.run();

	    // global clean-up...
        // fclaw_app_destroy will fail if there has been an abort
        if(!dirty_memory)
        {
            fclaw_app_destroy (app);
        }

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