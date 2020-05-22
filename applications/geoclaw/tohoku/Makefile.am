## Makefile.am in applications/geoclaw/tohoku

bin_PROGRAMS += applications/geoclaw/tohoku/tohoku

applications_geoclaw_tohoku_tohoku_SOURCES = \
	applications/geoclaw/tohoku/tohoku_user.cpp \
	applications/geoclaw/tohoku/tohoku_user.h \
	applications/geoclaw/tohoku/tohoku.cpp


applications_geoclaw_tohoku_tohoku_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_tohoku_tohoku_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
