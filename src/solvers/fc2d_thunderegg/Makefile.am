# -*- Makefile -*-
libthunderegg_generated_headers =

libthunderegg_installed_headers = \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg.h \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg_options.h \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg_physical_bc.h \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_starpatch.h \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_fivepoint.h \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_varpoisson.h \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_heat.h

libthunderegg_compiled_sources = \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg.cpp \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg_options.c \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg_vector.cpp \
	src/solvers/fc2d_thunderegg/fc2d_thunderegg_physical_bc.c \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_starpatch.cpp \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_fivepoint.cpp \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_varpoisson.cpp \
	src/solvers/fc2d_thunderegg/operators/fc2d_thunderegg_heat.cpp \
	src/solvers/fc2d_thunderegg/fort_4.6/thunderegg_apply_bc_default.f90 \
	src/solvers/fc2d_thunderegg/fort_4.6/thunderegg_eval_bc_default.f90

## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_thunderegg/libfc2d_thundergg.la

## Named variables that can be referenced from other libraries/apps
FCLAW_THUNDEREGG_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_thunderegg \
	                       -I@top_srcdir@/src/solvers/fc2d_thunderegg/operators \
						   $(THUNDEREGG_CFLAGS)

FCLAW_THUNDEREGG_LDADD    = @top_builddir@/src/solvers/fc2d_thunderegg/libfc2d_thundergg.la \
						   $(THUNDEREGG_LIBS)

## Sources needed to build this library
src_solvers_fc2d_thunderegg_libfc2d_thundergg_la_SOURCES = \
          $(libthunderegg_internal_headers) \
          $(libthunderegg_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_thunderegg_libfc2d_thundergg_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_THUNDEREGG_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_THUNDEREGG_CPPFLAGS)

LDADD        += $(FCLAW_THUNDEREGG_LDADD)

include_HEADERS += $(libthunderegg_installed_headers)
