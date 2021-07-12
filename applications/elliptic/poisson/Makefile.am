## Makefile.am in applications/elliptic/poisson

bin_PROGRAMS += applications/elliptic/poisson/poisson

applications_elliptic_poisson_poisson_SOURCES = \
	applications/elliptic/poisson/poisson.cpp \
	applications/elliptic/poisson/poisson_user.cpp \
	applications/elliptic/poisson/poisson_user.h \
	applications/elliptic/poisson/poisson_options.c \
	applications/elliptic/poisson/poisson_diagnostics.c \
	applications/elliptic/poisson/hsmooth.f90 \
	applications/elliptic/poisson/fortran/poisson_setprob.f90 \
	applications/elliptic/poisson/fortran/poisson_fort_rhs.f90 \
	applications/elliptic/poisson/fortran/poisson_compute_error.f \
	applications/elliptic/poisson/fortran/poisson_qexact.f90 \
	applications/elliptic/poisson/fortran/poisson_apply_bc.f90 \
	applications/elliptic/poisson/fortran/poisson_eval_bc.f90 \
	applications/elliptic/poisson/fortran/poisson_tag4refinement.f \
	applications/elliptic/poisson/fortran/poisson_tag4coarsening.f \
	applications/elliptic/poisson/fortran/poisson_output.f90


EXTRA_DIST += \
	applications/elliptic/poisson/fclaw_options.ini

#	applications/elliptic/poisson/regression.ini \
#	applications/elliptic/poisson/afterframe.m \
#	applications/elliptic/poisson/setplot2.m \
#	applications/elliptic/poisson/poisson.ipynb \
#	applications/elliptic/poisson/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_poisson_poisson_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_poisson_poisson_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
