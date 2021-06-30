## Makefile.am in applications/clawpack/shallow/2d/radialdam

bin_PROGRAMS += applications/clawpack/shallow/2d/radialdam/radialdam

applications_clawpack_shallow_2d_radialdam_radialdam_SOURCES = \
	applications/clawpack/shallow/2d/radialdam/radialdam.cpp \
	applications/clawpack/shallow/2d/radialdam/radialdam_options.c \
	applications/clawpack/shallow/2d/radialdam/radialdam_user.cpp \
	applications/clawpack/shallow/2d/radialdam/radialdam_user.h \
	applications/clawpack/shallow/2d/radialdam/fclaw2d_map_pillowdisk.c \
	applications/clawpack/shallow/2d/radialdam/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/shallow/2d/radialdam/setprob.f \
	applications/clawpack/shallow/2d/radialdam/gradient_exceeds_th.f90 \
	applications/clawpack/shallow/2d/radialdam/rpn2cons_update_manifold.f \
	applications/clawpack/shallow/2d/radialdam/fdisc.f \
	applications/clawpack/shallow/2d/radialdam/user_4.6/qinit.f \
	applications/clawpack/shallow/2d/radialdam/user_4.6/setaux.f90 \
	applications/clawpack/shallow/2d/radialdam/user_4.6/tag4refinement.f \
	applications/clawpack/shallow/2d/radialdam/user_4.6/tag4coarsening.f \
	applications/clawpack/shallow/2d/radialdam/user_5.0/qinit.f90 \
	applications/clawpack/shallow/2d/radialdam/user_5.0/setaux.f90 \
	applications/clawpack/shallow/2d/radialdam/user_5.0/tag4refinement.f90 \
	applications/clawpack/shallow/2d/radialdam/user_5.0/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/rp/shallow_user_fort.h \
	applications/clawpack/shallow/2d/rp/clawpack46_rpn2.f \
	applications/clawpack/shallow/2d/rp/clawpack46_rpt2.f \
	applications/clawpack/shallow/2d/rp/clawpack46_rpn2_manifold.f90 \
	applications/clawpack/shallow/2d/rp/clawpack46_rpt2_manifold.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpn2_manifold.f90 \
	applications/clawpack/shallow/2d/rp/clawpack5_rpt2_manifold.f90

EXTRA_DIST += \
	applications/clawpack/shallow/2d/radialdam/fclaw_options.ini \
	applications/clawpack/shallow/2d/radialdam/regression.ini \
	applications/clawpack/shallow/2d/radialdam/afterframe.m \
	applications/clawpack/shallow/2d/radialdam/setplot2.m \
	applications/clawpack/shallow/2d/radialdam/mapc2m.m \
	applications/clawpack/shallow/2d/radialdam/map1d.m

applications_clawpack_shallow_2d_radialdam_radialdam_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_shallow_2d_radialdam_radialdam_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
