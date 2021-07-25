## Makefile.am in applications/cudaclaw/euler/2d/shockbubble

bin_PROGRAMS += applications/cudaclaw/euler/2d/shockbubble_cuda/shockbubble


applications_cudaclaw_euler_2d_shockbubble_cuda_shockbubble_SOURCES = \
	applications/cudaclaw/euler/2d/shockbubble_cuda/shockbubble.cpp \
	applications/cudaclaw/euler/2d/shockbubble_cuda/shockbubble_options.c \
	applications/cudaclaw/euler/2d/shockbubble_cuda/shockbubble_user.cpp \
	applications/cudaclaw/euler/2d/shockbubble_cuda/shockbubble_user.h \
	applications/cudaclaw/euler/2d/shockbubble_cuda/fdisc.f \
	applications/cudaclaw/euler/2d/shockbubble_cuda/setprob.f90 \
	applications/cudaclaw/euler/2d/shockbubble_cuda/user_cuda/rp2euler.cu \
	applications/cudaclaw/euler/2d/shockbubble_cuda/user_4.6/qinit.f90 \
	applications/cudaclaw/euler/2d/shockbubble_cuda/user_4.6/bc2.f \
	applications/cudaclaw/euler/2d/shockbubble_cuda/user_4.6/src2.f

EXTRA_DIST += \
	applications/cudaclaw/euler/2d/shockbubble_cuda/fclaw_options.ini \
	applications/cudaclaw/euler/2d/shockbubble_cuda/regression.ini \
	applications/cudaclaw/euler/2d/shockbubble_cuda/afterframe.m \
	applications/cudaclaw/euler/2d/shockbubble_cuda/setplot2.m

applications_cudaclaw_euler_2d_shockbubble_cuda_shockbubble_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

SHOCKBUBBLE_DIR = applications/cudaclaw/euler/2d/shockbubble_cuda

applications_cudaclaw_euler_2d_shockbubble_cuda_shockbubble_LDADD = \
        $(LDADD) \
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



