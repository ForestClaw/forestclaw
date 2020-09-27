## Makefile.am in applications/geoclaw/chile2010

bin_PROGRAMS += applications/geoclaw/chile2010/chile2010

applications_geoclaw_chile2010_chile2010_SOURCES = \
	applications/geoclaw/chile2010/chile2010_user.h \
	applications/geoclaw/chile2010/chile2010.cpp


applications_geoclaw_chile2010_chile2010_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_GEOCLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_geoclaw_chile2010_chile2010_LDADD = \
         $(LDADD) \
         $(FCLAW_GEOCLAW_LDADD) \
         $(FCLAW_CLAWPATCH_LDADD) \
         $(FCLAW_LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
