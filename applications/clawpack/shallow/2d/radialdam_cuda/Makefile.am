## Makefile.am in applications/clawpack/shallow/2d/radialdam_cuda

bin_PROGRAMS += applications/clawpack/shallow/2d/radialdam_cuda/radialdam

applications_clawpack_shallow_2d_radialdam_cuda_radialdam_SOURCES = \
	applications/clawpack/shallow/2d/radialdam_cuda/radialdam.cpp \
	applications/clawpack/shallow/2d/radialdam_cuda/radialdam_options.c \
	applications/clawpack/shallow/2d/radialdam_cuda/radialdam_user.cpp \
	applications/clawpack/shallow/2d/radialdam_cuda/radialdam_user.h \
	applications/clawpack/shallow/2d/radialdam_cuda/fclaw2d_map_pillowdisk.c \
	applications/clawpack/shallow/2d/radialdam_cuda/setprob.f \
	applications/clawpack/shallow/2d/radialdam_cuda/fdisc.f \
	applications/clawpack/shallow/2d/radialdam_cuda/user_4.6/qinit.f \
	applications/clawpack/shallow/2d/radialdam_cuda/user_5.0/qinit.f90 \
	applications/clawpack/shallow/2d/radialdam_cuda/user_5.0/setaux.f90 \
	applications/clawpack/shallow/2d/radialdam_cuda/user_5.0/tag4refinement.f90 \
	applications/clawpack/shallow/2d/radialdam_cuda/user_5.0/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/rp/shallow_user_fort.h \
	applications/clawpack/shallow/2d/rp/clawpack46_rpn2.f \
	applications/clawpack/shallow/2d/rp/clawpack46_rpt2.f \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2_manifold.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2_manifold.f90 \
	applications/clawpack/shallow/2d/radialdam_cuda/user_cuda/rpn2shallow.cu

EXTRA_DIST += \
	applications/clawpack/shallow/2d/radialdam_cuda/fclaw_options.ini \
	applications/clawpack/shallow/2d/radialdam_cuda/regression.ini \
	applications/clawpack/shallow/2d/radialdam_cuda/afterframe.m \
	applications/clawpack/shallow/2d/radialdam_cuda/setplot2.m \
	applications/clawpack/shallow/2d/radialdam_cuda/mapc2m.m \
	applications/clawpack/shallow/2d/radialdam_cuda/map1d.m

applications_clawpack_shallow_2d_radialdam_cuda_radialdam_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

RADIALDAM_DIR = applications/clawpack/shallow/2d/radialdam_cuda

applications_clawpack_shallow_2d_radialdam_cuda_radialdam_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD) \
		$(RADIALDAM_DIR)/user_cuda/devicelink.o

## nvcc -dlink has to b called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
$(RADIALDAM_DIR)/user_cuda/devicelink.o : \
    $(RADIALDAM_DIR)/user_cuda/rpn2shallow.o \
    src/solvers/fc2d_cudaclaw/libcudaclaw.la 
		$(NVCC) -dlink  \
		$(RADIALDAM_DIR)/user_cuda/rpn2shallow.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(RADIALDAM_DIR)/user_cuda/devicelink.o 
