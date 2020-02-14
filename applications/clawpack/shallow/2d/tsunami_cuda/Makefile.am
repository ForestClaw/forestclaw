bin_PROGRAMS += applications/clawpack/shallow/2d/tsunami_cuda/tsunami

applications_clawpack_shallow_2d_tsunami_cuda_tsunami_SOURCES = \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami.cpp \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami_options.c \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami_user.cpp \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami_user.h \
	applications/clawpack/shallow/2d/tsunami_cuda/setprob.f \
	applications/clawpack/shallow/2d/tsunami_cuda/bathy.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/qinit.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/tag4refinement.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/setaux.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_cuda/rpn2_tsunami.cu \
	applications/clawpack/shallow/2d/rp/rpn2_tsunami.f90 \
	applications/clawpack/shallow/2d/rp/rpt2_tsunami.f90 \
	applications/clawpack/shallow/2d/rp/geoclaw_riemann_utils.f90

EXTRA_DIST += \
	applications/clawpack/shallow/2d/tsunami_cuda/fclaw_options.ini \
	applications/clawpack/shallow/2d/tsunami_cuda/regression.ini \
	applications/clawpack/shallow/2d/tsunami_cuda/afterframe.m \
	applications/clawpack/shallow/2d/tsunami_cuda/setplot2.m

TSUNAMI_DIR = applications/clawpack/shallow/2d/tsunami_cuda

applications_clawpack_shallow_2d_tsunami_cuda_tsunami_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_shallow_2d_tsunami_cuda_tsunami_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_LDADD) \
		$(TSUNAMI_DIR)/user_cuda/devicelink.o

## nvcc -dlink has to b called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
$(TSUNAMI_DIR)/user_cuda/devicelink.o : \
    $(TSUNAMI_DIR)/user_cuda/rpn2_tsunami.o \
    src/solvers/fc2d_cudaclaw/libcudaclaw.la 
		$(NVCC) -dlink  \
		$(TSUNAMI_DIR)/user_cuda/rpn2_tsunami.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(TSUNAMI_DIR)/user_cuda/devicelink.o 

