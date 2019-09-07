## Makefile.am in applications/elliptic/mgtest

bin_PROGRAMS += applications/elliptic/mgtest/mgtest

applications_elliptic_mgtest_mgtest_SOURCES = \
	applications/elliptic/mgtest/mgtest.cpp \
	applications/elliptic/mgtest/mgtest_user.cpp \
	applications/elliptic/mgtest/mgtest_user.h \
	applications/elliptic/mgtest/mgtest_options.c \
	applications/elliptic/mgtest/fortran/mgtest_setprob.f90 \
	applications/elliptic/mgtest/fortran/mgtest_fort_rhs.f90 \
	applications/elliptic/mgtest/fortran/mgtest_compute_error.f \
	applications/elliptic/mgtest/fortran/mgtest_qexact.f90 \
	applications/elliptic/mgtest/fortran/mgtest_apply_bc.f90 \
	applications/elliptic/mgtest/fortran/mgtest_eval_bc.f90 \
	applications/elliptic/mgtest/fortran/mgtest_tag4refinement.f90 \
	applications/elliptic/mgtest/fortran/mgtest_tag4coarsening.f90 \
	applications/elliptic/mgtest/fortran/mgtest_output.f90


EXTRA_DIST += \
	applications/elliptic/mgtest/fclaw_options.ini

#	applications/elliptic/mgtest/regression.ini \
#	applications/elliptic/mgtest/afterframe.m \
#	applications/elliptic/mgtest/setplot2.m \
#	applications/elliptic/mgtest/mgtest.ipynb \
#	applications/elliptic/mgtest/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

applications_elliptic_mgtest_mgtest_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_MULTIGRID_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_elliptic_mgtest_mgtest_LDADD = \
        $(LDADD) \
        $(FCLAW_MULTIGRID_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
