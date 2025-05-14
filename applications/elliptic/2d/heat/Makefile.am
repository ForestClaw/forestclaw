## Makefile.am in applications/elliptic/2d/heat

bin_PROGRAMS += applications/elliptic/2d/heat/heat

TESTS += applications/elliptic/2d/heat/regressions.sh

applications_elliptic_2d_heat_heat_SOURCES = \
	applications/elliptic/2d/heat/heat.cpp \
	applications/elliptic/2d/heat/heat_run.c \
	applications/elliptic/2d/heat/heat_user.cpp \
	applications/elliptic/2d/heat/heat_user.h \
	applications/elliptic/2d/heat/heat_options.c \
	applications/elliptic/2d/heat/heat_diagnostics.c \
	applications/elliptic/2d/heat/hsmooth.f90 \
	applications/elliptic/2d/heat/fortran/heat_setprob.f90 \
	applications/elliptic/2d/heat/fortran/heat_fort_rhs.f90 \
	applications/elliptic/2d/heat/fortran/heat_init.f90 \
	applications/elliptic/2d/heat/fortran/heat_compute_error.f \
	applications/elliptic/2d/heat/fortran/heat_qexact.f90 \
	applications/elliptic/2d/heat/fortran/heat_bc2.f \
	applications/elliptic/2d/heat/fortran/heat_apply_bc.f90 \
	applications/elliptic/2d/heat/fortran/heat_eval_bc.f90 \
	applications/elliptic/2d/heat/fortran/heat_tag4refinement.f \
	applications/elliptic/2d/heat/fortran/heat_tag4coarsening.f \
	applications/elliptic/2d/heat/fortran/heat_output.f90


EXTRA_DIST += \
	applications/elliptic/2d/heat/fclaw_options.ini

#	applications/elliptic/2d/heat/regression.ini \
#	applications/elliptic/2d/heat/afterframe.m \
#	applications/elliptic/2d/heat/setplot2.m \
#	applications/elliptic/2d/heat/heat.ipynb \
#	applications/elliptic/2d/heat/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_2d_heat_heat_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_2d_heat_heat_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
