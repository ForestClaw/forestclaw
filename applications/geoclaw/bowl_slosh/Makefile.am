## Makefile.am in applications/geoclaw/bowl_slosh

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/geoclaw/bowl_slosh/slosh

applications_geoclaw_bowl_slosh_slosh_SOURCES = \
	applications/geoclaw/bowl_slosh/slosh_user.cpp \
	applications/geoclaw/bowl_slosh/slosh_user.h \
	applications/geoclaw/bowl_slosh/slosh.cpp \
	applications/geoclaw/bowl_slosh/fclaw2d_map_nomap.c


applications_geoclaw_bowl_slosh_slosh_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_bowl_slosh_slosh_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
