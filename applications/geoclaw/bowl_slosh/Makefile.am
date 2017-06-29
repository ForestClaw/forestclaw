## Makefile.am in applications/geoclaw/bowl_slosh

bin_PROGRAMS += applications/geoclaw/bowl_slosh/slosh

applications_geoclaw_bowl_slosh_slosh_SOURCES = \
	applications/geoclaw/bowl_slosh/slosh_user.cpp \
	applications/geoclaw/bowl_slosh/slosh_user.h \
	applications/geoclaw/bowl_slosh/slosh.cpp \
	applications/geoclaw/bowl_slosh/qinit.f90


applications_geoclaw_bowl_slosh_slosh_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_bowl_slosh_slosh_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
