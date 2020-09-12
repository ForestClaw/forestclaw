## Makefile.am in applications/geoclaw/bowl_slosh

bin_PROGRAMS += applications/geoclaw/bowl_slosh/slosh

applications_geoclaw_bowl_slosh_slosh_SOURCES = \
	applications/geoclaw/bowl_slosh/slosh_user.cpp \
	applications/geoclaw/bowl_slosh/slosh_user.h \
	applications/geoclaw/bowl_slosh/slosh.cpp \
	applications/geoclaw/bowl_slosh/slosh_options.c \
	applications/geoclaw/bowl_slosh/qinit.f90

EXTRA_DIST += \
	applications/geoclaw/bowl_slosh/fclaw_options.ini \
	applications/geoclaw/bowl_slosh/make_data.py \
	applications/geoclaw/bowl_slosh/make_plots.py \
	applications/geoclaw/bowl_slosh/make_topo.py \
	applications/geoclaw/bowl_slosh/afterframe.m \
	applications/geoclaw/bowl_slosh/bathy.m \
	applications/geoclaw/bowl_slosh/setplot2.m \
	applications/geoclaw/bowl_slosh/wetcell.m \
	applications/geoclaw/bowl_slosh/plot_region.m

applications_geoclaw_bowl_slosh_slosh_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_GEOCLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_geoclaw_bowl_slosh_slosh_LDADD = \
         $(LDADD) \
         $(FCLAW_GEOCLAW_LDADD) \
         $(FCLAW_CLAWPATCH_LDADD) \
         $(FCLAW_LDADD)
