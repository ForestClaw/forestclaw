## Makefile.am in applications/geoclaw/bowl_slosh

bin_PROGRAMS += applications/geoclaw/bowl_slosh/slosh

applications_geoclaw_bowl_slosh_slosh_SOURCES = \
	applications/geoclaw/bowl_slosh/slosh_user.cpp \
	applications/geoclaw/bowl_slosh/slosh_user.h \
	applications/geoclaw/bowl_slosh/slosh.cpp \
	applications/geoclaw/bowl_slosh/qinit.f90

EXTRA_DIST += \
	applications/geoclaw/bowl_slosh/fclaw_options.ini \
	applications/geoclaw/bowl_slosh/make_data.py \
	applications/geoclaw/bowl_slosh/make_plots.py \
	applications/geoclaw/bowl_slosh/make_topo.py \
	applications/geoclaw/bowl_slosh/afterframe.py \
	applications/geoclaw/bowl_slosh/bathy.py \
	applications/geoclaw/bowl_slosh/setplot2.py \
	applications/geoclaw/bowl_slosh/wetcell.py \
	applications/geoclaw/bowl_slosh/plot_region.py

applications_geoclaw_bowl_slosh_slosh_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_bowl_slosh_slosh_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
