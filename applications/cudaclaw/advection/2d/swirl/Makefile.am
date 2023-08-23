## Makefile.am in applications/cudaclaw/advection/2d/swirl_cuda

bin_PROGRAMS += applications/cudaclaw/advection/2d/swirl/swirl_cuda

#TESTS += applications/cudaclaw/advection/2d/swirl/regressions.txt

applications_cudaclaw_advection_2d_swirl_swirl_cuda_SOURCES = \
	applications/cudaclaw/advection/2d/swirl/swirl_user.cpp \
	applications/cudaclaw/advection/2d/swirl/swirl_options.c \
	applications/cudaclaw/advection/2d/swirl/swirl.cpp \
	applications/cudaclaw/advection/2d/swirl/psi.f \
	applications/cudaclaw/advection/2d/swirl/user_cuda/qinit.f \
	applications/cudaclaw/advection/2d/swirl/user_cuda/rpn2adv.cu \
	applications/cudaclaw/advection/2d/swirl/user_cuda/rpt2adv.cu \
	applications/cudaclaw/advection/2d/swirl/user_cuda/b4step2.cu


EXTRA_DIST += \
	applications/cudaclaw/advection/2d/swirl/fclaw_options.ini \
	applications/cudaclaw/advection/2d/swirl/regression.ini \
	applications/cudaclaw/advection/2d/swirl/afterframe.m \
	applications/cudaclaw/advection/2d/swirl/setplot2.m

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_cudaclaw_advection_2d_swirl_swirl_cuda_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

SWIRL_DIR = applications/cudaclaw/advection/2d/swirl

applications_cudaclaw_advection_2d_swirl_swirl_cuda_LDADD = \
        $(LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD) \
		$(SWIRL_DIR)/user_cuda/devicelink.o 

## nvcc -dlink has to b called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
$(SWIRL_DIR)/user_cuda/devicelink.o : \
	$(SWIRL_DIR)/user_cuda/rpn2adv.o \
	$(SWIRL_DIR)/user_cuda/rpt2adv.o \
	$(SWIRL_DIR)/user_cuda/b4step2.o \
	src/solvers/fc2d_cudaclaw/libcudaclaw.la
		$(NVCC) $(CUDA_CFLAGS) -dlink  \
		$(SWIRL_DIR)/user_cuda/rpn2adv.o \
		$(SWIRL_DIR)/user_cuda/rpt2adv.o \
		$(SWIRL_DIR)/user_cuda/b4step2.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(SWIRL_DIR)/user_cuda/devicelink.o 


