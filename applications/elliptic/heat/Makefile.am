## Makefile.am in applications/elliptic/heat

bin_PROGRAMS += applications/elliptic/heat/heat

applications_elliptic_heat_heat_SOURCES = \
	applications/elliptic/heat/heat.cpp \
	applications/elliptic/heat/heat_run.c \
	applications/elliptic/heat/heat_user.cpp \
	applications/elliptic/heat/heat_user.h \
	applications/elliptic/heat/heat_options.c \
	applications/elliptic/heat/heat_diagnostics.c \
	applications/elliptic/heat/hsmooth.f90 \
	applications/elliptic/heat/fortran/heat_setprob.f90 \
	applications/elliptic/heat/fortran/heat_fort_rhs.f90 \
	applications/elliptic/heat/fortran/heat_init.f90 \
	applications/elliptic/heat/fortran/heat_compute_error.f \
	applications/elliptic/heat/fortran/heat_qexact.f90 \
	applications/elliptic/heat/fortran/heat_apply_bc.f90 \
	applications/elliptic/heat/fortran/heat_eval_bc.f90 \
	applications/elliptic/heat/fortran/heat_tag4refinement.f90 \
	applications/elliptic/heat/fortran/heat_tag4coarsening.f90 \
	applications/elliptic/heat/fortran/heat_output.f90


EXTRA_DIST += \
	applications/elliptic/heat/fclaw_options.ini

#	applications/elliptic/heat/regression.ini \
#	applications/elliptic/heat/afterframe.m \
#	applications/elliptic/heat/setplot2.m \
#	applications/elliptic/heat/heat.ipynb \
#	applications/elliptic/heat/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_heat_heat_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_heat_heat_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
