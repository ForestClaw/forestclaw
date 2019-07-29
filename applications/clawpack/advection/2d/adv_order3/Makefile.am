## Makefile.am in applications/clawpack/advection/2d/adv_order3

bin_PROGRAMS += applications/clawpack/advection/2d/adv_order3/periodic

applications_clawpack_advection_2d_adv_order3_periodic_SOURCES = \
	applications/clawpack/advection/2d/adv_order3/periodic_user.cpp \
	applications/clawpack/advection/2d/adv_order3/periodic_user.h \
	applications/clawpack/advection/2d/adv_order3/periodic_options.c \
	applications/clawpack/advection/2d/adv_order3/periodic.cpp \
	applications/clawpack/advection/2d/adv_order3/setprob.f \
    applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/adv_order3/user_4.6/qinit.f \
	applications/clawpack/advection/2d/adv_order3/user_4.6/rpn2.f \
	applications/clawpack/advection/2d/adv_order3/user_4.6/rpt2.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/adv_order3/fclaw_options.ini \
	applications/clawpack/advection/2d/adv_order3/afterframe.m \
	applications/clawpack/advection/2d/adv_order3/setplot2.m \
	applications/clawpack/advection/2d/adv_order3/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_clawpack_advection_2d_adv_order3_periodic_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_2d_adv_order3_periodic_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
