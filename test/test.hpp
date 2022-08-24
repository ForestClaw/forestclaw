#ifndef FCLAW_TEST_HPP
#define FCLAW_TEST_HPP
#include <fclaw2d_domain.h>
#include <fclaw2d_options.h>
#include <doctest.h>
#include <csetjmp>
fclaw2d_domain_t* create_test_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt);

void fclaw_test_expect_abort();

std::jmp_buf& fclaw_test_get_jump_buffer();

#define CHECK_SC_ABORTED(expr) \
{ \
    bool aborted = false; \
	fclaw_test_expect_abort(); \
	if(!setjmp(fclaw_test_get_jump_buffer())){ \
	    expr; \
	}else{ \
		aborted=true; \
	} \
	CHECK_UNARY(aborted); \
}

#endif