## Makefile.am in applications/elliptic/2d/allencahn

bin_PROGRAMS += applications/elliptic/2d/allencahn/allencahn

TESTS += applications/elliptic/2d/allencahn/regressions.sh

applications_elliptic_2d_allencahn_allencahn_SOURCES = \
	applications/elliptic/2d/allencahn/allencahn.cpp \
	applications/elliptic/2d/allencahn/allencahn_run.c \
	applications/elliptic/2d/allencahn/allencahn_user.cpp \
	applications/elliptic/2d/allencahn/allencahn_user.h \
	applications/elliptic/2d/allencahn/allencahn_options.c \
	applications/elliptic/2d/allencahn/hsmooth.f90 \
	applications/elliptic/2d/allencahn/fortran/allencahn_setprob.f90 \
	applications/elliptic/2d/allencahn/fortran/allencahn_bc2.f \
	applications/elliptic/2d/allencahn/fortran/allencahn_fort_rhs.f90 \
	applications/elliptic/2d/allencahn/fortran/allencahn_qinit.f \
	applications/elliptic/2d/allencahn/fortran/allencahn_qexact.f90 \
	applications/elliptic/2d/allencahn/fortran/allencahn_apply_bc.f90 \
	applications/elliptic/2d/allencahn/fortran/allencahn_eval_bc.f90 \
	applications/elliptic/2d/allencahn/fortran/allencahn_tag4refinement.f \
	applications/elliptic/2d/allencahn/fortran/allencahn_tag4coarsening.f


EXTRA_DIST += \
	applications/elliptic/2d/allencahn/fclaw_options.ini

#	applications/elliptic/2d/allencahn/regression.ini \
#	applications/elliptic/2d/allencahn/afterframe.m \
#	applications/elliptic/2d/allencahn/setplot2.m \
#	applications/elliptic/2d/allencahn/allencahn.ipynb \
#	applications/elliptic/2d/allencahn/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_2d_allencahn_allencahn_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_2d_allencahn_allencahn_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
