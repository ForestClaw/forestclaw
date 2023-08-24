bin_PROGRAMS += applications/clawpack/advection/2d/replicated/replicated

TESTS += applications/clawpack/advection/2d/replicated/regressions.sh

applications_clawpack_advection_2d_replicated_replicated_SOURCES = \
	applications/clawpack/advection/2d/replicated/replicated.cpp \
	applications/clawpack/advection/2d/replicated/replicated_user.h \
    applications/clawpack/advection/2d/replicated/replicated_options.c \
	applications/clawpack/advection/2d/replicated/replicated_user.cpp \
	applications/clawpack/advection/2d/replicated/setprob.f \
	applications/clawpack/advection/2d/replicated/psi.f \
	applications/clawpack/advection/2d/replicated/fdisc.f \
	applications/clawpack/advection/2d/replicated/user_4.6/qinit.f \
	applications/clawpack/advection/2d/replicated/user_5.0/qinit.f90 \
    applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/replicated/fclaw_options.ini \
	applications/clawpack/advection/2d/replicated/regression.ini \
	applications/clawpack/advection/2d/replicated/ex_00001.ini \
	applications/clawpack/advection/2d/replicated/ex_00004.ini \
	applications/clawpack/advection/2d/replicated/ex_00016.ini \
	applications/clawpack/advection/2d/replicated/p_00001.py \
	applications/clawpack/advection/2d/replicated/p_00004.py \
	applications/clawpack/advection/2d/replicated/p_00016.py \
	applications/clawpack/advection/2d/replicated/afterframe.m \
	applications/clawpack/advection/2d/replicated/setplot2.m \
	applications/clawpack/advection/2d/replicated/mapc2m.m

applications_clawpack_advection_2d_replicated_replicated_CPPFLAGS = \
        $(AM_CPPFLAGS)  \
        $(FCLAW_CLAWPACK46_CPPFLAGS)  \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_replicated_replicated_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

