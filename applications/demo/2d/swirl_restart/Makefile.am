## Makefile.am in applications/demo/2d/swirl_restart

bin_PROGRAMS += applications/demo/2d/swirl_restart/swirl

TESTS += applications/demo/2d/swirl_restart/regressions.sh

applications_demo_2d_swirl_restart_swirl_SOURCES = \
	applications/demo/2d/swirl_restart/swirl_user.cpp \
	applications/demo/2d/swirl_restart/swirl_user.h \
	applications/demo/2d/swirl_restart/swirl_ray.c \
	applications/demo/2d/swirl_restart/swirl_options.c \
	applications/demo/2d/swirl_restart/swirl.cpp \
	applications/demo/2d/swirl_restart/psi.f \
	applications/demo/2d/swirl_restart/setprob.f \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/demo/2d/swirl_restart/user_4.6/qinit.f \
	applications/demo/2d/swirl_restart/user_4.6/b4step2.f \
	applications/demo/2d/swirl_restart/user_5.0/qinit.f90 \
	applications/demo/2d/swirl_restart/user_5.0/b4step2.f90

EXTRA_DIST += \
	applications/demo/2d/swirl_restart/fclaw_options.ini \
	applications/demo/2d/swirl_restart/regression.ini \
	applications/demo/2d/swirl_restart/afterframe.m \
	applications/demo/2d/swirl_restart/setplot2.m \
	applications/demo/2d/swirl_restart/swirl.ipynb \
	applications/demo/2d/swirl_restart/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_demo_2d_swirl_restart_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_demo_2d_swirl_restart_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
