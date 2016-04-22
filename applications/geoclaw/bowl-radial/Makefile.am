## Makefile.am in applications/geoclaw/bowl-radial

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/geoclaw/bowl-radial/bowl

applications_geoclaw_bowl-radial_bowl_SOURCES = \
	applications/geoclaw/bowl-radial/bowl_user.cpp \
	applications/geoclaw/bowl-radial/bowl_user.h \
	applications/geoclaw/bowl-radial/bowl.cpp \
	applications/geoclaw/bowl-radial/fclaw2d_map_nomap.c \


applications_geoclaw_bowl-radial_bowl_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_bowl-radial_bowl_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
