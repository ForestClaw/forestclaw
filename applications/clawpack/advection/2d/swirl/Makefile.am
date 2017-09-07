## Makefile.am in applications/clawpack/advection/2d/swirl

bin_PROGRAMS += applications/clawpack/advection/2d/swirl/swirl

applications_clawpack_advection_2d_swirl_swirl_SOURCES = \
	applications/clawpack/advection/2d/swirl/swirl_user.cpp \
	applications/clawpack/advection/2d/swirl/swirl_user.h \
	applications/clawpack/advection/2d/swirl/swirl.cpp \
	applications/clawpack/advection/2d/swirl/psi.f \
	applications/clawpack/advection/2d/swirl/setprob.f \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/swirl/user_4.6/qinit.f \
	applications/clawpack/advection/2d/swirl/user_4.6/b4step2.f \
	applications/clawpack/advection/2d/swirl/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/swirl/user_5.0/b4step2.f90

applications_clawpack_advection_2d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_clawpack_advection_2d_swirl_swirl_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
