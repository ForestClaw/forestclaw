## Makefile.am in applications/cudaclaw/shallow/2d/bump

bin_PROGRAMS += applications/cudaclaw/shallow/2d/bump/bump

applications_cudaclaw_shallow_2d_bump_bump_SOURCES = \
	applications/cudaclaw/shallow/2d/bump/bump.cpp \
	applications/cudaclaw/shallow/2d/bump/bump_options.c \
	applications/cudaclaw/shallow/2d/bump/bump_user.cpp \
	applications/cudaclaw/shallow/2d/bump/bump_user.h \
	applications/cudaclaw/shallow/2d/bump/user_cuda/qinit.f \
	applications/cudaclaw/shallow/2d/bump/user_cuda/rpn2shallow.cu

EXTRA_DIST += \
	applications/cudaclaw/shallow/2d/bump/fclaw_options.ini \
	applications/cudaclaw/shallow/2d/bump/regression.ini \
	applications/cudaclaw/shallow/2d/bump/afterframe.m \
	applications/cudaclaw/shallow/2d/bump/setplot2.m \
	applications/cudaclaw/shallow/2d/bump/mapc2m.m

applications_cudaclaw_shallow_2d_bump_bump_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

BUMP_DIR = applications/cudaclaw/shallow/2d/bump

applications_cudaclaw_shallow_2d_bump_bump_LDADD = \
        $(LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD) \
		$(BUMP_DIR)/user_cuda/devicelink.o

## nvcc -dlink has to b called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
$(BUMP_DIR)/user_cuda/devicelink.o : \
    $(BUMP_DIR)/user_cuda/rpn2shallow.o \
    src/solvers/fc2d_cudaclaw/libcudaclaw.la 
		$(NVCC) -dlink  \
		$(BUMP_DIR)/user_cuda/rpn2shallow.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(BUMP_DIR)/user_cuda/devicelink.o 


