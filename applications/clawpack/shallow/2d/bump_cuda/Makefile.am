## Makefile.am in applications/clawpack/shallow/2d/bump_cuda

bin_PROGRAMS += applications/clawpack/shallow/2d/bump_cuda/bump

applications_clawpack_shallow_2d_bump_cuda_bump_SOURCES = \
	applications/clawpack/shallow/2d/bump_cuda/bump.cpp \
	applications/clawpack/shallow/2d/bump_cuda/bump_options.c \
	applications/clawpack/shallow/2d/bump_cuda/bump_user.cpp \
	applications/clawpack/shallow/2d/bump_cuda/bump_user.h \
	applications/clawpack/shallow/2d/bump_cuda/fclaw2d_map_fivepatch.c \
	applications/clawpack/shallow/2d/bump_cuda/setprob.f \
	applications/clawpack/shallow/2d/bump_cuda/user_4.6/qinit.f \
	applications/clawpack/shallow/2d/bump_cuda/user_4.6/tag4refinement.f90 \
	applications/clawpack/shallow/2d/bump_cuda/user_4.6/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/bump_cuda/user_5.0/qinit.f90 \
	applications/clawpack/shallow/2d/bump_cuda/user_5.0/setaux.f90 \
	applications/clawpack/shallow/2d/bump_cuda/user_5.0/tag4refinement.f90 \
	applications/clawpack/shallow/2d/bump_cuda/user_5.0/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/rp/shallow_user_fort.h \
	applications/clawpack/shallow/2d/rp/clawpack46_rpn2.f \
	applications/clawpack/shallow/2d/rp/clawpack46_rpt2.f \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2_manifold.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2_manifold.f90 \
	applications/clawpack/shallow/2d/radialdam_cuda/user_cuda/rpn2shallow.cu

EXTRA_DIST += \
	applications/clawpack/shallow/2d/bump_cuda/fclaw_options.ini \
	applications/clawpack/shallow/2d/bump_cuda/regression.ini \
	applications/clawpack/shallow/2d/bump_cuda/afterframe.m \
	applications/clawpack/shallow/2d/bump_cuda/setplot2.m \
	applications/clawpack/shallow/2d/bump_cuda/mapc2m.m

applications_clawpack_shallow_2d_bump__cuda_bump_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

BUMP_DIR = applications/clawpack/shallow/2d/bump_cuda

applications_clawpack_shallow_2d_bump_cuda_bump_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD) \
		$(BUMPDIR_DIR)/user_cuda/devicelink.o
