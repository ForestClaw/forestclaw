## Makefile.am in applications/clawpack/shallow/2d/bump

bin_PROGRAMS += applications/clawpack/shallow/2d/bump/bump

applications_clawpack_shallow_2d_bump_bump_SOURCES = \
	applications/clawpack/shallow/2d/bump/bump.cpp \
	applications/clawpack/shallow/2d/bump/bump_options.c \
	applications/clawpack/shallow/2d/bump/bump_user.cpp \
	applications/clawpack/shallow/2d/bump/bump_user.h \
	applications/clawpack/shallow/2d/bump/fclaw2d_map_pillowdisk.c \
	applications/clawpack/shallow/2d/bump/setprob.f \
	applications/clawpack/shallow/2d/bump/fdisc.f \
	applications/clawpack/shallow/2d/bump/user_4.6/qinit.f \
	applications/clawpack/shallow/2d/bump/user_4.6/tag4refinement.f90 \
	applications/clawpack/shallow/2d/bump/user_4.6/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/bump/user_5.0/qinit.f90 \
	applications/clawpack/shallow/2d/bump/user_5.0/setaux.f90 \
	applications/clawpack/shallow/2d/bump/user_5.0/tag4refinement.f90 \
	applications/clawpack/shallow/2d/bump/user_5.0/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/rp/shallow_user_fort.h \
	applications/clawpack/shallow/2d/rp/clawpack46_rpn2.f \
	applications/clawpack/shallow/2d/rp/clawpack46_rpt2.f \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2_manifold.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2_manifold.f90

EXTRA_DIST += \
	applications/clawpack/shallow/2d/bump/fclaw_options.ini \
	applications/clawpack/shallow/2d/bump/regression.ini \
	applications/clawpack/shallow/2d/bump/afterframe.m \
	applications/clawpack/shallow/2d/bump/setplot2.m \
	applications/clawpack/shallow/2d/bump/mapc2m.m \
	applications/clawpack/shallow/2d/bump/map1d.m

applications_clawpack_shallow_2d_bump_bump_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_shallow_2d_bump_bump_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
