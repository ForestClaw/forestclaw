## Makefile.am in applications/clawpack/euler/2d/shockbubble

bin_PROGRAMS += applications/clawpack/euler/2d/shockbubble/shockbubble

TESTS += applications/clawpack/euler/2d/shockbubble/regressions.sh

applications_clawpack_euler_2d_shockbubble_shockbubble_SOURCES = \
	applications/clawpack/euler/2d/shockbubble/shockbubble.cpp \
	applications/clawpack/euler/2d/shockbubble/shockbubble_options.c \
	applications/clawpack/euler/2d/shockbubble/shockbubble_user.cpp \
	applications/clawpack/euler/2d/shockbubble/shockbubble_user.h \
	applications/clawpack/euler/2d/shockbubble/fdisc.f \
	applications/clawpack/euler/2d/shockbubble/setprob.f \
	applications/clawpack/euler/2d/shockbubble/difference_exceeds_th.f90 \
	applications/clawpack/euler/2d/shockbubble/user_4.6/qinit.f \
	applications/clawpack/euler/2d/shockbubble/user_4.6/setaux.f \
	applications/clawpack/euler/2d/shockbubble/user_4.6/bc2.f \
	applications/clawpack/euler/2d/shockbubble/user_4.6/src2.f \
	applications/clawpack/euler/2d/shockbubble/user_4.6/tag4coarsening.f \
	applications/clawpack/euler/2d/shockbubble/user_4.6/tag4refinement.f \
	applications/clawpack/euler/2d/shockbubble/user_5.0/qinit.f90 \
	applications/clawpack/euler/2d/shockbubble/user_5.0/setaux.f90 \
	applications/clawpack/euler/2d/shockbubble/user_5.0/bc2.f90 \
	applications/clawpack/euler/2d/shockbubble/user_5.0/src2.f90 \
	applications/clawpack/euler/2d/shockbubble/user_5.0/tag4coarsening.f \
	applications/clawpack/euler/2d/shockbubble/user_5.0/tag4refinement.f \
	applications/clawpack/euler/2d/rp/euler_user_fort.h \
	applications/clawpack/euler/2d/rp/clawpack46_rpn2_euler4.f \
	applications/clawpack/euler/2d/rp/clawpack46_rpt2_euler4.f \
	applications/clawpack/euler/2d/rp/clawpack46_rpn2_euler5.f \
	applications/clawpack/euler/2d/rp/clawpack46_rpt2_euler5.f \
	applications/clawpack/euler/2d/rp/clawpack5_rpn2_euler4.f90 \
	applications/clawpack/euler/2d/rp/clawpack5_rpt2_euler4.f90 \
	applications/clawpack/euler/2d/rp/clawpack5_rpn2_euler5.f90 \
	applications/clawpack/euler/2d/rp/clawpack5_rpt2_euler5.f90

EXTRA_DIST += \
	applications/clawpack/euler/2d/shockbubble/fclaw_options.ini \
	applications/clawpack/euler/2d/shockbubble/regression.ini \
	applications/clawpack/euler/2d/shockbubble/afterframe.m \
	applications/clawpack/euler/2d/shockbubble/setplot2.m

applications_clawpack_euler_2d_shockbubble_shockbubble_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_euler_2d_shockbubble_shockbubble_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)


