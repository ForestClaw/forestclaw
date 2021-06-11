## Makefile.am in applications/clawpack/advection/2d/periodic

bin_PROGRAMS += applications/clawpack/advection/2d/periodic/periodic

applications_clawpack_advection_2d_periodic_periodic_SOURCES = \
	applications/clawpack/advection/2d/periodic/periodic_user.cpp \
	applications/clawpack/advection/2d/periodic/periodic_user.h \
	applications/clawpack/advection/2d/periodic/periodic_options.c \
	applications/clawpack/advection/2d/periodic/periodic.cpp \
	applications/clawpack/advection/2d/periodic/setprob.f \
	applications/clawpack/advection/2d/periodic/fdisc.f \
    applications/clawpack/advection/2d/periodic/user_4.6/clawpack46_setaux.f \
	applications/clawpack/advection/2d/periodic/user_4.6/qinit.f \
	applications/clawpack/advection/2d/periodic/user_4.6/tag4refinement.f \
	applications/clawpack/advection/2d/periodic/user_4.6/tag4coarsening.f \
    applications/clawpack/advection/2d/periodic/user_5.0/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/periodic/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/periodic/user_5.0/tag4refinement.f \
	applications/clawpack/advection/2d/periodic/user_5.0/tag4coarsening.f \
    applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/periodic/fclaw_options.ini \
	applications/clawpack/advection/2d/periodic/regression.ini \
	applications/clawpack/advection/2d/periodic/afterframe.m \
	applications/clawpack/advection/2d/periodic/setplot2.m \
	applications/clawpack/advection/2d/periodic/periodic.ipynb \
	applications/clawpack/advection/2d/periodic/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_clawpack_advection_2d_periodic_periodic_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_2d_periodic_periodic_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
