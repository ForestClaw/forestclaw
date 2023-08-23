## Makefile.am in applications/elliptic/allencahn

bin_PROGRAMS += applications/elliptic/allencahn/allencahn

TESTS += applications/elliptic/allencahn/regressions.txt

applications_elliptic_allencahn_allencahn_SOURCES = \
	applications/elliptic/allencahn/allencahn.cpp \
	applications/elliptic/allencahn/allencahn_run.c \
	applications/elliptic/allencahn/allencahn_user.cpp \
	applications/elliptic/allencahn/allencahn_user.h \
	applications/elliptic/allencahn/allencahn_options.c \
	applications/elliptic/allencahn/hsmooth.f90 \
	applications/elliptic/allencahn/fortran/allencahn_setprob.f90 \
	applications/elliptic/allencahn/fortran/allencahn_bc2.f \
	applications/elliptic/allencahn/fortran/allencahn_fort_rhs.f90 \
	applications/elliptic/allencahn/fortran/allencahn_qinit.f \
	applications/elliptic/allencahn/fortran/allencahn_qexact.f90 \
	applications/elliptic/allencahn/fortran/allencahn_apply_bc.f90 \
	applications/elliptic/allencahn/fortran/allencahn_eval_bc.f90 \
	applications/elliptic/allencahn/fortran/allencahn_tag4refinement.f \
	applications/elliptic/allencahn/fortran/allencahn_tag4coarsening.f


EXTRA_DIST += \
	applications/elliptic/allencahn/fclaw_options.ini

#	applications/elliptic/allencahn/regression.ini \
#	applications/elliptic/allencahn/afterframe.m \
#	applications/elliptic/allencahn/setplot2.m \
#	applications/elliptic/allencahn/allencahn.ipynb \
#	applications/elliptic/allencahn/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_allencahn_allencahn_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_allencahn_allencahn_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
