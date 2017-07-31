## Makefile.am in applications/geoclaw/teton

bin_PROGRAMS += applications/geoclaw/teton/teton

applications_geoclaw_teton_teton_SOURCES = \
	applications/geoclaw/teton/teton_user.cpp \
	applications/geoclaw/teton/teton_user.h \
	applications/geoclaw/teton/qinit.f90 \
	applications/geoclaw/teton/teton.cpp


applications_geoclaw_teton_teton_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_teton_teton_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
