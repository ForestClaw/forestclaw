## Makefile.am in applications/cudaclaw/advection/2d/swirl

bin_PROGRAMS += applications/cudaclaw/advection/2d/swirl/swirl

applications_cudaclaw_advection_2d_swirl_swirl_SOURCES = \
	applications/cudaclaw/advection/2d/swirl/swirl_user.cpp \
	applications/cudaclaw/advection/2d/swirl/swirl_user.h \
	applications/cudaclaw/advection/2d/swirl/swirl_options.c \
	applications/cudaclaw/advection/2d/swirl/swirl.cpp \
	applications/cudaclaw/advection/2d/swirl/psi.f \
	applications/cudaclaw/advection/2d/swirl/setprob.f \
	applications/cudaclaw/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/cudaclaw/advection/2d/swirl/user_5.0/qinit.f90 \
	applications/cudaclaw/advection/2d/swirl/user_5.0/b4step2.f90

EXTRA_DIST += \
	applications/cudaclaw/advection/2d/swirl/fclaw_options.ini \
	applications/cudaclaw/advection/2d/swirl/regression.ini \
	applications/cudaclaw/advection/2d/swirl/afterframe.m \
	applications/cudaclaw/advection/2d/swirl/setplot2.m \
	applications/cudaclaw/advection/2d/swirl/swirl.ipynb \
	applications/cudaclaw/advection/2d/swirl/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_cudaclaw_advection_2d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CUDACLAW5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_cudaclaw_advection_2d_swirl_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CUDACLAW5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
