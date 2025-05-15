# -*- Makefile -*-
libthunderegg_3d_generated_headers =

libthunderegg_3d_installed_headers = \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg.h \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_options.h \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_physical_bc.h

libthunderegg_3d_compiled_sources = \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg.cpp \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_options.c \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_vector.cpp \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_physical_bc.c \
	src/solvers/fc3d_thunderegg/fort_4.6/thunderegg_apply_bc_default.f90 \
	src/solvers/fc3d_thunderegg/fort_4.6/thunderegg_eval_bc_default.f90

## Name of library to build
lib_LTLIBRARIES += src/solvers/fc3d_thunderegg/libfc3d_thunderegg.la

## Named variables that can be referenced from other libraries/apps
FCLAW_THUNDEREGG_3D_CPPFLAGS = -I@top_srcdir@/src/solvers/fc3d_thunderegg \
						       $(THUNDEREGG_CFLAGS)

FCLAW_THUNDEREGG_3D_LDADD    = @top_builddir@/src/solvers/fc3d_thunderegg/libfc3d_thunderegg.la \
						   	   $(THUNDEREGG_LIBS)

## Sources needed to build this library
src_solvers_fc3d_thunderegg_libfc3d_thunderegg_la_SOURCES = \
          $(libthunderegg_3d_internal_headers) \
          $(libthunderegg_3d_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are
## included in AM_CPPFLAGS, but are included here for clarity (although makes
## for lengthy compile lines)
src_solvers_fc3d_thunderegg_libfc3d_thunderegg_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_THUNDEREGG_3D_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_3D_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_THUNDEREGG_3D_CPPFLAGS)

LDADD        += $(FCLAW_THUNDEREGG_3D_LDADD)

include_HEADERS += $(libthunderegg_3d_installed_headers)

## UNIT TESTS
check_PROGRAMS += src/solvers/fc3d_thunderegg/fc3d_thunderegg.TEST
TESTS += src/solvers/fc3d_thunderegg/fc3d_thunderegg.TEST

src_solvers_fc3d_thunderegg_fc3d_thunderegg_TEST_SOURCES = \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg.h.TEST.cpp \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_options.h.TEST.cpp \
	src/solvers/fc3d_thunderegg/fc3d_thunderegg_vector_TEST.cpp

src_solvers_fc3d_thunderegg_fc3d_thunderegg_TEST_CPPFLAGS = \
    $(test_libtestutils_la_CPPFLAGS) \
    $(src_solvers_fc3d_thunderegg_libfc3d_thunderegg_la_CPPFLAGS)

src_solvers_fc3d_thunderegg_fc3d_thunderegg_TEST_LDADD = \
	test/libtestutils.la \
	$(test_libtestutils_la_LDADD) \
    $(LDADD) \
    $(FCLAW_THUNDEREGG_3D_LDADD) \
    $(FCLAW_CLAWPATCH_LDADD) \
    $(FCLAW_LDADD)
