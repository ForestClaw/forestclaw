bin_PROGRAMS += applications/clawpack/shallow/2d/tsunami_cuda/tsunami

applications_clawpack_shallow_2d_tsunami_cuda_tsunami_SOURCES = \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami.cpp \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami_options.c \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami_user.cpp \
	applications/clawpack/shallow/2d/tsunami_cuda/tsunami_user.h \
	applications/clawpack/shallow/2d/tsunami_cuda/setprob.f \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/qinit.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/tag4refinement.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/bathy.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/user_4.6/setaux.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/rp/rpn2_tsunami.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/rp/rpt2_tsunami.f90 \
	applications/clawpack/shallow/2d/tsunami_cuda/rp/rpn2_tsunami.cu \
	applications/clawpack/shallow/2d/tsunami_cuda/rp/geoclaw_riemann_utils.f90

EXTRA_DIST += \
	applications/clawpack/shallow/2d/tsunami_cuda/fclaw_options.ini \
	applications/clawpack/shallow/2d/tsunami_cuda/regression.ini \
	applications/clawpack/shallow/2d/tsunami_cuda/afterframe.m \
	applications/clawpack/shallow/2d/tsunami_cuda/setplot2.m

applications_clawpack_shallow_2d_tsunami_cuda_tsunami_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_shallow_2d_tsunami_cuda_tsunami_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
