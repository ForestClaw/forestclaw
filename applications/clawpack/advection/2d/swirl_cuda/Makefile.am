## Makefile.am in applications/clawpack/advection/2d/swirl_cuda

bin_PROGRAMS += applications/clawpack/advection/2d/swirl_cuda/swirl

applications_clawpack_advection_2d_swirl_cuda_swirl_SOURCES = \
	applications/clawpack/advection/2d/swirl_cuda/swirl_user.cpp \
	applications/clawpack/advection/2d/swirl_cuda/swirl_user.h \
	applications/clawpack/advection/2d/swirl_cuda/swirl_options.c \
	applications/clawpack/advection/2d/swirl_cuda/swirl.cpp \
	applications/clawpack/advection/2d/swirl_cuda/psi.f \
	applications/clawpack/advection/2d/swirl_cuda/setprob.f \
    applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/swirl_cuda/user_4.6/qinit.f \
	applications/clawpack/advection/2d/swirl_cuda/user_4.6/b4step2.f \
	applications/clawpack/advection/2d/swirl_cuda/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/swirl_cuda/user_5.0/b4step2.f90 \
	applications/clawpack/advection/2d/swirl_cuda/user_cuda/qinit.f \
	applications/clawpack/advection/2d/swirl_cuda/user_cuda/rpn2adv.cu \
	applications/clawpack/advection/2d/swirl_cuda/user_cuda/rpt2adv.cu \
	applications/clawpack/advection/2d/swirl_cuda/user_cuda/b4step2.cu


EXTRA_DIST += \
	applications/clawpack/advection/2d/swirl_cuda/fclaw_options.ini \
	applications/clawpack/advection/2d/swirl_cuda/regression.ini \
	applications/clawpack/advection/2d/swirl_cuda/afterframe.m \
	applications/clawpack/advection/2d/swirl_cuda/setplot2.m \
	applications/clawpack/advection/2d/swirl_cuda/swirl.ipynb \
	applications/clawpack/advection/2d/swirl_cuda/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_clawpack_advection_2d_swirl_cuda_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

SWIRL_DIR = applications/clawpack/advection/2d/swirl_cuda

applications_clawpack_advection_2d_swirl_cuda_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
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
	src/solvers/fc2d_cudaclaw/libcudaclaw.la \
		$(NVCC) $(CUDA_CFLAGS) -dlink  \
		$(SWIRL_DIR)/user_cuda/rpn2adv.o \
		$(SWIRL_DIR)/user_cuda/rpt2adv.o \
		$(SWIRL_DIR)/user_cuda/b4step2.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(SWIRL_DIR)/user_cuda/devicelink.o 


