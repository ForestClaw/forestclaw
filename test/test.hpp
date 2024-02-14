#ifndef FCLAW_TEST_HPP
#define FCLAW_TEST_HPP
#include <doctest.h>
#include <csetjmp>
#include <fclaw_base.h>

bool test_output_vtk();
bool test_indirect();

void fclaw_test_expect_abort();
void fclaw_test_clear_expect_abort();

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
	fclaw_test_clear_expect_abort(); \
}

#endif

extern "C"
{
#define TEST_2D_CLAWPATCH46_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(test_2d_clawpatch46_fort_interpolate_face, \
                             TEST_2D_CLAWPATCH46_FORT_INTERPOLATE_FACE)
void TEST_2D_CLAWPATCH46_FORT_INTERPOLATE_FACE(const int* mx, const int* my, 
                                               const int* mbc,const int* meqn,
                                               double qcoarse[],double qfine[],
                                               const int* idir, const int* iside,
                                               const int* num_neighbors,
                                               const int* refratio, const int* igrid,
                                               struct fclaw_patch_transform_data** 
                                               transform_cptr);

#define TEST_2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(test_2d_clawpatch46_fort_interpolate_corner, \
                     TEST_2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER)
void TEST_2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, 
                                                 const int* mbc,const int* meqn, 
                                                 const int* a_refratio, double this_q[],
                                                 double neighbor_q[], const int* a_corner,
                                                 struct fclaw_patch_transform_data** 
                                                 transform_cptr);
}