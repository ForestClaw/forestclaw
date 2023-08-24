## Makefile.am in applications/geoclaw/bowl_radial_slosh

bin_PROGRAMS += applications/geoclaw/bowl_radial_slosh/bowl_radial_slosh

TESTS += applications/geoclaw/bowl_radial_slosh/regressions.sh

applications_geoclaw_bowl_radial_slosh_bowl_radial_slosh_SOURCES = \
	applications/geoclaw/bowl_radial_slosh/slosh/slosh_user.cpp \
	applications/geoclaw/bowl_radial_slosh/slosh/slosh_user.h \
	applications/geoclaw/bowl_radial_slosh/slosh/slosh_options.c \
	applications/geoclaw/bowl_radial_slosh/slosh/qinit.f90 \
	applications/geoclaw/bowl_radial_slosh/radial/radial_user.cpp \
	applications/geoclaw/bowl_radial_slosh/radial/radial_user.h \
	applications/geoclaw/bowl_radial_slosh/radial/radial_options.c \
	applications/geoclaw/bowl_radial_slosh/bowl.cpp

EXTRA_DIST += \
	applications/geoclaw/bowl_raidal_slosh/raidal/fclaw_options.ini \
	applications/geoclaw/bowl_raidal_slosh/raidal/make_data.py \
	applications/geoclaw/bowl_raidal_slosh/raidal/make_plots.py \
	applications/geoclaw/bowl_raidal_slosh/raidal/make_topo.py \
	applications/geoclaw/bowl_raidal_slosh/raidal/afterframe.m \
	applications/geoclaw/bowl_raidal_slosh/raidal/setplot2.m \
	applications/geoclaw/bowl_raidal_slosh/raidal/bathy.m \
	applications/geoclaw/bowl_raidal_slosh/raidal/bathy1d.m \
	applications/geoclaw/bowl_raidal_slosh/raidal/map1d.m \
	applications/geoclaw/bowl_raidal_slosh/raidal/compare_gauges.m \
	applications/geoclaw/bowl_radial_slosh/slosh/fclaw_options.ini \
	applications/geoclaw/bowl_radial_slosh/slosh/make_data.py \
	applications/geoclaw/bowl_radial_slosh/slosh/make_plots.py \
	applications/geoclaw/bowl_radial_slosh/slosh/make_topo.py \
	applications/geoclaw/bowl_radial_slosh/slosh/afterframe.m \
	applications/geoclaw/bowl_radial_slosh/slosh/bathy.m \
	applications/geoclaw/bowl_radial_slosh/slosh/setplot2.m \
	applications/geoclaw/bowl_radial_slosh/slosh/wetcell.m \
	applications/geoclaw/bowl_radial_slosh/slosh/plot_region.m

applications_geoclaw_bowl_radial_slosh_bowl_radial_slosh_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_GEOCLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_geoclaw_bowl_radial_slosh_bowl_radial_slosh_LDADD = \
         $(LDADD) \
         $(FCLAW_GEOCLAW_LDADD) \
         $(FCLAW_CLAWPATCH_LDADD) \
         $(FCLAW_LDADD)
