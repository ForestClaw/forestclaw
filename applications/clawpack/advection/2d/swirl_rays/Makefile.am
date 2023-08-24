## Makefile.am in applications/clawpack/advection/2d/swirl_rays

bin_PROGRAMS += applications/clawpack/advection/2d/swirl_rays/swirl

TESTS += applications/clawpack/advection/2d/swirl_rays/regressions.sh

applications_clawpack_advection_2d_swirl_rays_swirl_SOURCES = \
	applications/clawpack/advection/2d/swirl_rays/swirl_user.cpp \
	applications/clawpack/advection/2d/swirl_rays/swirl_user.h \
	applications/clawpack/advection/2d/swirl_rays/swirl_options.c \
	applications/clawpack/advection/2d/swirl_rays/swirl.cpp \
	applications/clawpack/advection/2d/swirl_rays/psi.f \
	applications/clawpack/advection/2d/swirl_rays/setprob.f \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/swirl_rays/user_4.6/qinit.f \
	applications/clawpack/advection/2d/swirl_rays/user_4.6/b4step2.f \
	applications/clawpack/advection/2d/swirl_rays/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/swirl_rays/user_5.0/b4step2.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/swirl_rays/fclaw_options.ini \
	applications/clawpack/advection/2d/swirl_rays/regression.ini \
	applications/clawpack/advection/2d/swirl_rays/afterframe.m \
	applications/clawpack/advection/2d/swirl_rays/setplot2.m \
	applications/clawpack/advection/2d/swirl_rays/swirl.ipynb \
	applications/clawpack/advection/2d/swirl_rays/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_clawpack_advection_2d_swirl_rays_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_2d_swirl_rays_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
