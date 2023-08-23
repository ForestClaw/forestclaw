## Makefile.am in applications/elliptic/heat_phasefield

bin_PROGRAMS += applications/elliptic/heat_phasefield/heat_phasefield

TESTS += applications/elliptic/heat_phasefield/regressions.txt

applications_elliptic_heat_phasefield_heat_phasefield_SOURCES = \
	applications/elliptic/heat_phasefield/heat_phasefield.cpp \
	applications/elliptic/heat_phasefield/heat/heat.cpp \
	applications/elliptic/heat_phasefield/heat/heat_run.c \
	applications/elliptic/heat_phasefield/heat/heat_user.cpp \
	applications/elliptic/heat_phasefield/heat/heat_user.h \
	applications/elliptic/heat_phasefield/heat/heat_options.c \
	applications/elliptic/heat_phasefield/heat/heat_diagnostics.c \
	applications/elliptic/heat_phasefield/heat/hsmooth.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_setprob.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_fort_rhs.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_init.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_compute_error.f \
	applications/elliptic/heat_phasefield/heat/fortran/heat_qexact.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_bc2.f \
	applications/elliptic/heat_phasefield/heat/fortran/heat_apply_bc.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_eval_bc.f90 \
	applications/elliptic/heat_phasefield/heat/fortran/heat_tag4refinement.f \
	applications/elliptic/heat_phasefield/heat/fortran/heat_tag4coarsening.f \
	applications/elliptic/heat_phasefield/heat/fortran/heat_output.f90 \
	applications/elliptic/heat_phasefield/phasefield/phasefield.cpp \
	applications/elliptic/heat_phasefield/phasefield/phasefield_run.c \
	applications/elliptic/heat_phasefield/phasefield/phasefield_operator.cpp \
	applications/elliptic/heat_phasefield/phasefield/phasefield_operator.h \
	applications/elliptic/heat_phasefield/phasefield/phasefield_patch_operator.cpp \
	applications/elliptic/heat_phasefield/phasefield/phasefield_patch_operator.h \
	applications/elliptic/heat_phasefield/phasefield/phasefield_user.cpp \
	applications/elliptic/heat_phasefield/phasefield/phasefield_user.h \
	applications/elliptic/heat_phasefield/phasefield/phasefield_options.c \
	applications/elliptic/heat_phasefield/phasefield/cellave.f \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_setprob.f90 \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_bc2.f \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_fort_rhs.f90 \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_qinit.f90 \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_apply_bc.f90 \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_eval_bc.f90 \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_tag4refinement.f \
	applications/elliptic/heat_phasefield/phasefield/fortran/phasefield_tag4coarsening.f




EXTRA_DIST += \
	applications/elliptic/heat_phasefield/fclaw_options.ini

#	applications/elliptic/heat_phasefield/regression.ini \
#	applications/elliptic/heat_phasefield/afterframe.m \
#	applications/elliptic/heat_phasefield/setplot2.m \
#	applications/elliptic/heat_phasefield/heat_phasefield.ipynb \
#	applications/elliptic/heat_phasefield/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_heat_phasefield_heat_phasefield_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_heat_phasefield_heat_phasefield_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
