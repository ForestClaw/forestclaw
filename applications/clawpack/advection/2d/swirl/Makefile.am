## Makefile.am in applications/clawpack/advection/2d/swirl

bin_PROGRAMS += applications/clawpack/advection/2d/swirl/swirl

applications_clawpack_advection_2d_swirl_swirl_SOURCES = \
	applications/clawpack/advection/2d/swirl/swirl_user.cpp \
	applications/clawpack/advection/2d/swirl/swirl_user.h \
	applications/clawpack/advection/2d/swirl/swirl_options.c \
	applications/clawpack/advection/2d/swirl/swirl.cpp \
	applications/clawpack/advection/2d/swirl/psi.f \
	applications/clawpack/advection/2d/swirl/setprob.f \
    applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/swirl/user_4.6/qinit.f \
	applications/clawpack/advection/2d/swirl/user_4.6/b4step2.f \
	applications/clawpack/advection/2d/swirl/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/swirl/user_5.0/b4step2.f90 \
	applications/clawpack/advection/2d/swirl/user_cuda/qinit.f \
	applications/clawpack/advection/2d/swirl/user_cuda/rpn2adv.cu \
	applications/clawpack/advection/2d/swirl/user_cuda/rpt2adv.cu \
	applications/clawpack/advection/2d/swirl/user_cuda/b4step2.cu


EXTRA_DIST += \
	applications/clawpack/advection/2d/swirl/fclaw_options.ini \
	applications/clawpack/advection/2d/swirl/regression.ini \
	applications/clawpack/advection/2d/swirl/afterframe.m \
	applications/clawpack/advection/2d/swirl/setplot2.m \
	applications/clawpack/advection/2d/swirl/swirl.ipynb \
	applications/clawpack/advection/2d/swirl/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_clawpack_advection_2d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_2d_swirl_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD) \
		applications/clawpack/advection/2d/swirl/user_cuda/devicelink.o 

## nvcc -dlink has to b called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
applications/clawpack/advection/2d/swirl/user_cuda/devicelink.o : \
	applications/clawpack/advection/2d/swirl/user_cuda/rpn2adv.o \
	applications/clawpack/advection/2d/swirl/user_cuda/rpt2adv.o \
	applications/clawpack/advection/2d/swirl/user_cuda/b4step2.o \
	src/solvers/fc2d_cudaclaw/libcudaclaw.la 
		$(NVCC) $(CUDA_CFLAGS) -dlink  \
		applications/clawpack/advection/2d/swirl/user_cuda/rpn2adv.o \
		applications/clawpack/advection/2d/swirl/user_cuda/rpt2adv.o \
		applications/clawpack/advection/2d/swirl/user_cuda/b4step2.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o applications/clawpack/advection/2d/swirl/user_cuda/devicelink.o 


