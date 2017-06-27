## Makefile.am in applications/geoclaw/bowl_radial

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/geoclaw/bowl_radial/bowl

applications_geoclaw_bowl_radial_bowl_SOURCES = \
	applications/geoclaw/bowl_radial/bowl_user.cpp \
	applications/geoclaw/bowl_radial/bowl_user.h \
	applications/geoclaw/bowl_radial/bowl.cpp


applications_geoclaw_bowl_radial_bowl_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_bowl_radial_bowl_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
