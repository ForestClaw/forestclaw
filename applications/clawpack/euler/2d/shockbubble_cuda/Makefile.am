## Makefile.am in applications/clawpack/euler/2d/shockbubble

bin_PROGRAMS += applications/clawpack/euler/2d/shockbubble_cuda/shockbubble


applications_clawpack_euler_2d_shockbubble_cuda_shockbubble_SOURCES = \
	applications/clawpack/euler/2d/shockbubble_cuda/shockbubble.cpp \
	applications/clawpack/euler/2d/shockbubble_cuda/shockbubble_options.c \
	applications/clawpack/euler/2d/shockbubble_cuda/shockbubble_user.cpp \
	applications/clawpack/euler/2d/shockbubble_cuda/shockbubble_user.h \
	applications/clawpack/euler/2d/shockbubble_cuda/fdisc.f \
	applications/clawpack/euler/2d/shockbubble_cuda/setprob.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_cuda/rp2euler.cu \
	applications/clawpack/euler/2d/shockbubble_cuda/user_4.6/qinit.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_4.6/setaux.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_4.6/bc2.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_4.6/src2.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_4.6/tag4coarsening.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_4.6/tag4refinement.f \
	applications/clawpack/euler/2d/shockbubble_cuda/user_5.0/qinit.f90 \
	applications/clawpack/euler/2d/shockbubble_cuda/user_5.0/setaux.f90 \
	applications/clawpack/euler/2d/shockbubble_cuda/user_5.0/bc2.f90 \
	applications/clawpack/euler/2d/shockbubble_cuda/user_5.0/src2.f90 \
	applications/clawpack/euler/2d/shockbubble_cuda/user_5.0/tag4coarsening.f90 \
	applications/clawpack/euler/2d/shockbubble_cuda/user_5.0/tag4refinement.f90 \
	applications/clawpack/euler/2d/rp/euler_user_fort.h \
	applications/clawpack/euler/2d/rp/clawpack46_rpn2_euler4.f \
	applications/clawpack/euler/2d/rp/clawpack46_rpt2_euler4.f \
	applications/clawpack/euler/2d/rp/clawpack46_rpn2_euler5.f \
	applications/clawpack/euler/2d/rp/clawpack46_rpt2_euler5.f \
	applications/clawpack/euler/2d/rp/clawpack5_rpn2_euler4.f90 \
	applications/clawpack/euler/2d/rp/clawpack5_rpt2_euler4.f90 \
	applications/clawpack/euler/2d/rp/clawpack5_rpn2_euler5.f90 \
	applications/clawpack/euler/2d/rp/clawpack5_rpt2_euler5.f90

EXTRA_DIST += \
	applications/clawpack/euler/2d/shockbubble_cuda/fclaw_options.ini \
	applications/clawpack/euler/2d/shockbubble_cuda/regression.ini \
	applications/clawpack/euler/2d/shockbubble_cuda/afterframe.m \
	applications/clawpack/euler/2d/shockbubble_cuda/setplot2.m

applications_clawpack_euler_2d_shockbubble_cuda_shockbubble_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

SHOCKBUBBLE_DIR = applications/clawpack/euler/2d/shockbubble_cuda

applications_clawpack_euler_2d_shockbubble_cuda_shockbubble_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_LDADD) \
		$(SHOCKBUBBLE_DIR)/user_cuda/devicelink.o


## nvcc -dlink has to b called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
$(SHOCKBUBBLE_DIR)/user_cuda/devicelink.o : \
    $(SHOCKBUBBLE_DIR)/user_cuda/rp2euler.o \
    src/solvers/fc2d_cudaclaw/libcudaclaw.la 
		$(NVCC) -dlink  \
		$(SHOCKBUBBLE_DIR)/user_cuda/rp2euler.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(SHOCKBUBBLE_DIR)/user_cuda/devicelink.o 



