## Makefile.am in applications/geoclaw/flat_bathymetry

bin_PROGRAMS += applications/geoclaw/flat_bathymetry/flatbathy

applications_geoclaw_flat_bathymetry_flatbathy_SOURCES = \
	applications/geoclaw/flat_bathymetry/flatbathy_user.cpp \
	applications/geoclaw/flat_bathymetry/flatbathy_user.h \
	applications/geoclaw/flat_bathymetry/flatbathy.cpp \
	applications/geoclaw/flat_bathymetry/qinit.f


applications_geoclaw_flat_bathymetry_flatbathy_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

applications_geoclaw_flat_bathymetry_flatbathy_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la
