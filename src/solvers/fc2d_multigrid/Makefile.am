# -*- Makefile -*-
libmultigrid_generated_headers =

libmultigrid_installed_headers = \
	src/solvers/fc2d_multigrid/fc2d_multigrid.h \
	src/solvers/fc2d_multigrid/fc2d_multigrid_options.h

libmultigrid_compiled_sources = \
	src/solvers/fc2d_multigrid/fc2d_multigrid.cpp \
	src/solvers/fc2d_multigrid/fc2d_multigrid_options.c \
	src/solvers/fc2d_multigrid/fc2d_multigrid_solve.c


## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_multigrid/libmultigrid.la

## Named variables that can be referenced from other libraries/apps
FCLAW_MULTIGRID_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_multigrid 
FCLAW_MULTIGRID_LDADD    = @top_builddir@/src/solvers/fc2d_multigrid/libmultigrid.la

## Sources needed to build this library
src_solvers_fc2d_multigrid_libmultigrid_la_SOURCES = \
          $(libmultigrid_internal_headers) \
          $(libmultigrid_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_multigrid_libmultigrid_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_MULTIGRID_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_MULTIGRID_CPPFLAGS)

LDADD        += $(FCLAW_MULTIGRID_LDADD)

include_HEADERS += $(libmultigrid_installed_headers)
