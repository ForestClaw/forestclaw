## Makefile.am in applications/elliptic/2d/phasefield

bin_PROGRAMS += applications/elliptic/2d/phasefield/phasefield

TESTS += applications/elliptic/2d/phasefield/regressions.sh

applications_elliptic_2d_phasefield_phasefield_SOURCES = \
	applications/elliptic/2d/phasefield/phasefield.cpp \
	applications/elliptic/2d/phasefield/phasefield_run.c \
	applications/elliptic/2d/phasefield/phasefield_operator.cpp \
	applications/elliptic/2d/phasefield/phasefield_operator.h \
	applications/elliptic/2d/phasefield/phasefield_patch_operator.cpp \
	applications/elliptic/2d/phasefield/phasefield_patch_operator.h \
	applications/elliptic/2d/phasefield/phasefield_user.cpp \
	applications/elliptic/2d/phasefield/phasefield_user.h \
	applications/elliptic/2d/phasefield/phasefield_options.c \
	applications/elliptic/2d/phasefield/cellave.f \
	applications/elliptic/2d/phasefield/fortran/phasefield_setprob.f90 \
	applications/elliptic/2d/phasefield/fortran/phasefield_bc2.f \
	applications/elliptic/2d/phasefield/fortran/phasefield_fort_rhs.f90 \
	applications/elliptic/2d/phasefield/fortran/phasefield_qinit.f90 \
	applications/elliptic/2d/phasefield/fortran/phasefield_apply_bc.f90 \
	applications/elliptic/2d/phasefield/fortran/phasefield_eval_bc.f90 \
	applications/elliptic/2d/phasefield/fortran/phasefield_tag4refinement.f \
	applications/elliptic/2d/phasefield/fortran/phasefield_tag4coarsening.f



EXTRA_DIST += \
	applications/elliptic/2d/phasefield/fclaw_options.ini

#	applications/elliptic/2d/phasefield/regression.ini \
#	applications/elliptic/2d/phasefield/afterframe.m \
#	applications/elliptic/2d/phasefield/setplot2.m \
#	applications/elliptic/2d/phasefield/phasefield.ipynb \
#	applications/elliptic/2d/phasefield/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_2d_phasefield_phasefield_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_elliptic_2d_phasefield_phasefield_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)