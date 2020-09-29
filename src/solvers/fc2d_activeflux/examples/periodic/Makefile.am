## Makefile.am in src/solvers/fc2d_activeflux

bin_PROGRAMS += src/solvers/fc2d_activeflux/examples/periodic/periodic

src_solvers_fc2d_activeflux_examples_periodic_periodic_SOURCES = \
	src/solvers/fc2d_activeflux/examples/periodic/periodic_user.cpp \
	src/solvers/fc2d_activeflux/examples/periodic/periodic_user.h \
	src/solvers/fc2d_activeflux/examples/periodic/periodic_options.c \
	src/solvers/fc2d_activeflux/examples/periodic/periodic.cpp \
	src/solvers/fc2d_activeflux/examples/periodic/setprob.f \
	src/solvers/fc2d_activeflux/examples/periodic/user_4.6/qinit.f \
	src/solvers/fc2d_activeflux/examples/periodic/user_4.6/compute_error.f \
	src/solvers/fc2d_activeflux/examples/periodic/user_4.6/tag4refinement.f \
	src/solvers/fc2d_activeflux/examples/periodic/user_4.6/tag4coarsening.f \
	src/solvers/fc2d_activeflux/examples/rp/rpn2_adv.f \
	src/solvers/fc2d_activeflux/examples/rp/rpt2_adv.f \
	src/solvers/fc2d_activeflux/examples/rp/rpn2cons_update_adv.f \
    applications/clawpack/advection/2d/all/advection_user_fort.h

EXTRA_DIST += \
	src/solvers/fc2d_activeflux/examples/periodic/fclaw_options.ini \
	src/solvers/fc2d_activeflux/examples/periodic/afterframe.m \
	src/solvers/fc2d_activeflux/examples/periodic/setplot2.m \
	src/solvers/fc2d_activeflux/examples/periodic/make_plots.py

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
src_solvers_fc2d_activeflux_examples_periodic_periodic_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_ACTIVEFLUX_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


src_solvers_fc2d_activeflux_examples_periodic_periodic_LDADD = \
        $(LDADD) \
        $(FCLAW_ACTIVEFLUX_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
