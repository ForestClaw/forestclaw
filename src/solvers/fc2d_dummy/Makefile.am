# -*- Makefile -*-
libdummy_generated_headers =
libdummy_installed_headers = \
	src/solvers/fc2d_dummy/fc2d_dummy.H

libdummy_compiled_sources = \
	src/solvers/fc2d_dummy/fc2d_dummy.cpp

lib_LTLIBRARIES += src/solvers/fc2d_dummy/libdummy.la

src_solvers_fc2d_dummy_libdummy_la_SOURCES = \
  $(libdummy_internal_headers) \
  $(libdummy_compiled_sources)

src_solvers_fc2d_dummy_libdummy_la_CPPFLAGS = \
  $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

FCLAW_DUMMY_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_dummy
FCLAW_DUMMY_LDADD = @top_builddir@/src/solvers/fc2d_dummy/libdummy.la

include_HEADERS += $(libdummy_installed_headers)
