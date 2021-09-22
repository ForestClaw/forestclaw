## Makefile.am in applications/cudaclaw/euler/2d/shockbubble

bin_PROGRAMS += applications/cudaclaw/euler/2d/shockbubble/shockbubble


applications_cudaclaw_euler_2d_shockbubble_shockbubble_SOURCES = \
	applications/cudaclaw/euler/2d/shockbubble/shockbubble.cpp \
	applications/cudaclaw/euler/2d/shockbubble/shockbubble_options.c \
	applications/cudaclaw/euler/2d/shockbubble/shockbubble_user.cpp \
	applications/cudaclaw/euler/2d/shockbubble/shockbubble_user.h \
	applications/cudaclaw/euler/2d/shockbubble/fdisc.f \
	applications/cudaclaw/euler/2d/shockbubble/setprob.f90 \
	applications/cudaclaw/euler/2d/shockbubble/user_cuda/rp2euler.cu \
	applications/cudaclaw/euler/2d/shockbubble/user_cuda/qinit.f90 \
	applications/cudaclaw/euler/2d/shockbubble/user_cuda/bc2.f \
	applications/cudaclaw/euler/2d/shockbubble/user_cuda/src2.f

EXTRA_DIST += \
	applications/cudaclaw/euler/2d/shockbubble/fclaw_options.ini \
	applications/cudaclaw/euler/2d/shockbubble/regression.ini \
	applications/cudaclaw/euler/2d/shockbubble/afterframe.m \
	applications/cudaclaw/euler/2d/shockbubble/setplot2.m

applications_cudaclaw_euler_2d_shockbubble_shockbubble_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

SHOCKBUBBLE_DIR = applications/cudaclaw/euler/2d/shockbubble

applications_cudaclaw_euler_2d_shockbubble_shockbubble_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAWCLAW_LDADD) \
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



